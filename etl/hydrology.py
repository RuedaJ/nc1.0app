from __future__ import annotations
from pathlib import Path
import numpy as np
import xarray as xr
import geopandas as gpd
import rioxarray as rxr
import rasterio as rio
from rasterio.enums import Resampling
from shapely.geometry import box

# Shapely 2.x function (preferred) with 1.x fallback
try:
    from shapely import union_all as _union_all
except Exception:
    from shapely.ops import unary_union as _union_all  # function form


def _ensure_2d(da: xr.DataArray) -> xr.DataArray:
    """Squeeze a 3D [band,y,x] raster to 2D [y,x] if needed."""
    if da.ndim == 3:
        da = da.squeeze(drop=True)
    return da


def clip_to_aoi(raster_path: str | Path, aoi_geojson: str | Path) -> xr.DataArray:
    """Clip raster to AOI; error if there is no overlap. Keeps CRS/transform."""
    da = _ensure_2d(rxr.open_rasterio(raster_path, masked=True))
    aoi = gpd.read_file(aoi_geojson)
    if aoi.empty:
        raise ValueError("AOI is empty.")

    # Reproject AOI to raster CRS if necessary
    if aoi.crs is not None and da.rio.crs is not None and str(aoi.crs) != str(da.rio.crs):
        aoi_r = aoi.to_crs(da.rio.crs)
    else:
        aoi_r = aoi

    # Overlap guard
    rb = da.rio.bounds()
    inter = box(*rb).intersection(_union_all(list(aoi_r.geometry)))
    if inter.is_empty:
        raise ValueError(
            f"No overlap between raster and AOI. "
            f"Raster bounds={rb}, AOI bounds={tuple(aoi_r.total_bounds)}"
        )

    clipped = da.rio.clip(aoi_r.geometry, aoi_r.crs, drop=True)

    # Ensure nodata handled and spatial metadata preserved
    clipped = clipped.where(np.isfinite(clipped), np.nan)
    if getattr(clipped.rio, "crs", None) is None and getattr(da.rio, "crs", None) is not None:
        clipped = clipped.rio.write_crs(da.rio.crs)
    try:
        _ = clipped.rio.transform(recalc=True)
    except Exception:
        clipped = clipped.rio.write_transform(da.rio.transform(recalc=True))
    return clipped


def slope_from_dem(dem_path: str | Path) -> xr.DataArray:
    """Compute slope (degrees) from DEM using finite differences; preserves CRS/transform."""
    dem = _ensure_2d(rxr.open_rasterio(dem_path, masked=True)).where(lambda x: np.isfinite(x), np.nan)

    # resolution returns (xres, yres); use absolute spacing
    resx, resy = dem.rio.resolution()
    arr = dem.values.astype("float32")

    with np.errstate(invalid="ignore", divide="ignore"):
        gy, gx = np.gradient(arr, abs(resy), abs(resx))  # dZ/dy, dZ/dx
        slope_rad = np.arctan(np.hypot(gx, gy))
        slope_deg = np.degrees(slope_rad)

    slope = xr.DataArray(
        slope_deg,
        coords=dem.coords,
        dims=dem.dims,
        attrs=dem.attrs,
        name="slope_deg",
    )
    if dem.rio.crs is not None:
        slope = slope.rio.write_crs(dem.rio.crs, inplace=False)
        slope = slope.rio.write_transform(dem.rio.transform(recalc=True), inplace=False)
    return slope


def reproject_match(src: xr.DataArray, template: xr.DataArray,
                    resampling: Resampling = Resampling.nearest) -> xr.DataArray:
    """Reproject src to exactly match template grid/CRS/transform."""
    src = _ensure_2d(src)
    template = _ensure_2d(template)
    if src.rio.crs is None:
        raise ValueError("Source raster has no CRS; cannot reproject_match.")
    if template.rio.crs is None:
        raise ValueError("Template raster has no CRS; cannot reproject_match.")
    out = src.rio.reproject_match(template, resampling=resampling)
    return out.where(np.isfinite(out), np.nan)


def _norm_0_1(da: xr.DataArray, invert: bool = False) -> xr.DataArray:
    """Min-max normalize to [0,1]; optional inversion (1-x). Keeps CRS/transform."""
    arr = da.values.astype("float32")
    finite = np.isfinite(arr)
    if finite.any():
        vmin = float(np.nanmin(arr[finite]))
        vmax = float(np.nanmax(arr[finite]))
        if vmax > vmin:
            norm = (arr - vmin) / (vmax - vmin)
        else:
            norm = np.zeros_like(arr, dtype="float32")
    else:
        norm = np.zeros_like(arr, dtype="float32")
    if invert:
        norm = 1.0 - norm

    out = xr.DataArray(norm, coords=da.coords, dims=da.dims, attrs=da.attrs)
    if da.rio.crs is not None:
        out = out.rio.write_crs(da.rio.crs, inplace=False)
        out = out.rio.write_transform(da.rio.transform(recalc=True), inplace=False)
    return out


def infiltration_score(awc_da: xr.DataArray,
                       slope_da: xr.DataArray,
                       lc_da: xr.DataArray | None,
                       w_awc: float, w_slope: float, w_lc: float) -> xr.DataArray:
    """
    Compute infiltration score in [0,1]:
    - Higher AWC -> higher score
    - Flatter slope -> higher score (slope inverted)
    - Land cover optional (defaults to 0.5)
    """
    awc_n = _norm_0_1(awc_da, invert=False)
    slope_n = _norm_0_1(slope_da, invert=True)  # flatter is better
    if lc_da is None:
        lc_n = xr.full_like(awc_n, 0.5)
    else:
        lc_n = _norm_0_1(lc_da, invert=False)
        # align LC to AWC grid if needed
        if (lc_n.sizes != awc_n.sizes) or (lc_n.rio.transform() != awc_n.rio.transform()):
            lc_n = reproject_match(lc_n, awc_n, resampling=Resampling.bilinear)

    weights = np.array([w_awc, w_slope, w_lc], dtype="float32")
    sw = float(weights.sum()) if float(weights.sum()) > 0 else 1.0
    score_vals = (w_awc * awc_n.values + w_slope * slope_n.values + w_lc * lc_n.values) / sw
    score_vals = np.clip(score_vals, 0.0, 1.0).astype("float32")

    score = xr.DataArray(
        score_vals,
        coords=awc_n.coords,
        dims=awc_n.dims,
        attrs=awc_n.attrs,
        name="InfiltrationScore",
    )
    if awc_n.rio.crs is not None:
        score = score.rio.write_crs(awc_n.rio.crs, inplace=False)
        score = score.rio.write_transform(awc_n.rio.transform(recalc=True), inplace=False)
    return score
