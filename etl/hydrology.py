from pathlib import Path
import numpy as np
import xarray as xr
import rioxarray as rxr
import rasterio as rio
import geopandas as gpd
from shapely.geometry import box
from rioxarray.exceptions import NoDataInBounds
from rasterio.enums import Resampling

# ───────────────────── Utilities ─────────────────────
def _normalize01(da: xr.DataArray, method: str = "percentile", pmin=2, pmax=98) -> xr.DataArray:
    """Normalize to 0..1 ignoring NaNs; robust to outliers."""
    vals = da.values.astype("float64")
    mask = np.isfinite(vals)
    if not np.any(mask):
        return xr.full_like(da, np.nan)
    if method == "percentile":
        a = np.nanpercentile(vals, pmin)
        b = np.nanpercentile(vals, pmax)
    else:
        a = np.nanmin(vals)
        b = np.nanmax(vals)
    if not np.isfinite(a) or not np.isfinite(b) or a == b:
        return xr.zeros_like(da)
    out = (da - a) / (b - a)
    return out.clip(0, 1)

def _ensure_like(ref: xr.DataArray, da_or_vals, name: str = None) -> xr.DataArray:
    """Ensure output shares ref's dims/coords/CRS/transform."""
    if isinstance(da_or_vals, xr.DataArray):
        out = da_or_vals
        need_rewrap = (
            tuple(out.dims) != tuple(ref.dims)
            or any(d not in out.coords for d in ref.dims)
        )
        if need_rewrap:
            out = xr.DataArray(
                np.asarray(out),
                coords={ref.dims[0]: ref.coords[ref.dims[0]], ref.dims[1]: ref.coords[ref.dims[1]]},
                dims=ref.dims,
                attrs=ref.attrs,
                name=name or out.name,
            )
    else:
        out = xr.DataArray(
            np.asarray(da_or_vals),
            coords={ref.dims[0]: ref.coords[ref.dims[0]], ref.dims[1]: ref.coords[ref.dims[1]]},
            dims=ref.dims,
            attrs=ref.attrs,
            name=name,
        )
    if getattr(out.rio, "crs", None) is None and getattr(ref.rio, "crs", None) is not None:
        out = out.rio.write_crs(ref.rio.crs)
    try:
        _ = out.rio.transform(recalc=True)
    except Exception:
        out = out.rio.write_transform(ref.rio.transform(recalc=True))
    return out

# ───────────────────── Core functions ─────────────────────
def clip_to_aoi(raster_path: Path, aoi_geojson: Path) -> xr.DataArray:
    da = rxr.open_rasterio(raster_path, masked=True).squeeze()  # (y,x)
    aoi = gpd.read_file(aoi_geojson)

    # CRS checks
    if da.rio.crs is None:
        raise ValueError(f"Raster has no CRS: {raster_path}")
    if aoi.crs is None:
        raise ValueError(f"AOI has no CRS: {aoi_geojson}")
    if str(aoi.crs) != str(da.rio.crs):
        aoi = aoi.to_crs(da.rio.crs)

    # Pre-check intersection to avoid NoDataInBounds
    rb = da.rio.bounds()  # (minx, miny, maxx, maxy)
    raster_poly = box(*rb)
    aoi_union = aoi.geometry.unary_union
    if not raster_poly.intersects(aoi_union):
        raise ValueError(f"No overlap between raster and AOI. Raster bounds={rb}, AOI bounds={aoi.total_bounds}")

    # Clip with drop, fallback to keep shape
    try:
        clipped = da.rio.clip(aoi.geometry, aoi.crs, drop=True, invert=False)
    except NoDataInBounds:
        clipped = da.rio.clip(aoi.geometry, aoi.crs, drop=False, invert=False)
    return clipped

def slope_from_dem(dem_path: Path) -> xr.DataArray:
    """Compute slope in degrees from a DEM GeoTIFF/VRT using central differences."""
    dem = rxr.open_rasterio(dem_path, masked=True).squeeze()  # (y,x)
    # pixel size
    transform = dem.rio.transform()
    px = transform.a
    py = -transform.e  # transform.e is negative for north-up rasters
    dx = abs(px) if px else 1.0
    dy = abs(py) if py else 1.0

    z = dem.values.astype("float64")
    # handle all-nan
    if np.all(~np.isfinite(z)):
        return xr.full_like(dem, np.nan)

    # gradients
    dz_dy, dz_dx = np.gradient(z, dy, dx)  # rows->y, cols->x
    slope_rad = np.arctan(np.sqrt(dz_dx**2 + dz_dy**2))
    slope_deg = np.degrees(slope_rad)

    slope = xr.DataArray(
        slope_deg,
        coords={dem.dims[0]: dem.coords[dem.dims[0]], dem.dims[1]: dem.coords[dem.dims[1]]},
        dims=dem.dims,
        attrs=dem.attrs,
        name="Slope",
    )
    if getattr(slope.rio, "crs", None) is None:
        slope = slope.rio.write_crs(dem.rio.crs)
    try:
        _ = slope.rio.transform(recalc=True)
    except Exception:
        slope = slope.rio.write_transform(dem.rio.transform(recalc=True))
    return slope

def reproject_match(src_da: xr.DataArray, match_da: xr.DataArray, resampling=Resampling.bilinear) -> xr.DataArray:
    out = src_da.rio.reproject_match(match_da, resampling=resampling)
    return _ensure_like(match_da, out, name=src_da.name or "reprojected"

    )

def infiltration_score(awc_da: xr.DataArray, slope_da: xr.DataArray, lc_da: xr.DataArray | None,
                       w_awc: float, w_slp: float, w_lc: float) -> xr.DataArray:
    """Combine AWC, slope, and optional land cover into a 0..100 score."""
    ref = awc_da

    # Normalize inputs (robust percentiles)
    awc_n = _normalize01(awc_da, method="percentile")

    # Slope: low slope → higher infiltration; clamp at 30°
    slope_clamped = slope_da.clip(0, 30)
    slope_n = 1.0 - (slope_clamped / 30.0)  # 1 at 0°, 0 at 30°+
    slope_n = _ensure_like(ref, slope_n, name="SlopeIndex")

    # Land cover index
    if lc_da is None:
        lc_n = xr.full_like(ref, 0.5)
    else:
        lc_n = _normalize01(reproject_match(lc_da, ref, resampling=Resampling.nearest))

    # Weights normalized to sum 1
    w = np.array([w_awc, w_slp, w_lc], dtype="float64")
    if not np.isfinite(w).all() or w.sum() <= 0:
        w = np.array([1, 1, 1], dtype="float64")
    w = w / w.sum()

    # Weighted sum to 0..1, then scale to 0..100
    comp = (w[0] * _ensure_like(ref, awc_n) +
            w[1] * _ensure_like(ref, slope_n) +
            w[2] * _ensure_like(ref, lc_n))

    score01 = comp.clip(0, 1)
    score = (score01 * 100.0).astype("float32")
    out = _ensure_like(ref, score, name="InfiltrationScore")
    return out
