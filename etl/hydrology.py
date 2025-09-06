from pathlib import Path
import numpy as np
import xarray as xr
import rioxarray as rxr
import rasterio as rio
import geopandas as gpd
from shapely.geometry import box
from rioxarray.exceptions import NoDataInBounds
from rasterio.enums import Resampling

def _normalize01(da: xr.DataArray, method: str = "percentile", pmin=2, pmax=98) -> xr.DataArray:
    vals = da.values.astype("float64")
    mask = np.isfinite(vals)
    if not np.any(mask):
        return xr.full_like(da, np.nan)
    if method == "percentile":
        a = np.nanpercentile(vals, pmin)
        b = np.nanpercentile(vals, pmax)
    else:
        a = np.nanmin(vals); b = np.nanmax(vals)
    if not np.isfinite(a) or not np.isfinite(b) or a == b:
        return xr.zeros_like(da)
    return ((da - a) / (b - a)).clip(0, 1)

def _ensure_like(ref: xr.DataArray, da_or_vals, name: str | None = None) -> xr.DataArray:
    if isinstance(da_or_vals, xr.DataArray):
        out = da_or_vals
        need_rewrap = (tuple(out.dims) != tuple(ref.dims)) or any(d not in out.coords for d in ref.dims)
        if need_rewrap:
            out = xr.DataArray(np.asarray(out),
                               coords={ref.dims[0]: ref.coords[ref.dims[0]], ref.dims[1]: ref.coords[ref.dims[1]]},
                               dims=ref.dims, attrs=ref.attrs, name=name or out.name)
    else:
        out = xr.DataArray(np.asarray(da_or_vals),
                           coords={ref.dims[0]: ref.coords[ref.dims[0]], ref.dims[1]: ref.coords[ref.dims[1]]},
                           dims=ref.dims, attrs=ref.attrs, name=name)
    if getattr(out.rio, "crs", None) is None and getattr(ref.rio, "crs", None) is not None:
        out = out.rio.write_crs(ref.rio.crs)
    try:
        _ = out.rio.transform(recalc=True)
    except Exception:
        out = out.rio.write_transform(ref.rio.transform(recalc=True))
    return out

def clip_to_aoi(raster_path: Path, aoi_geojson: Path) -> xr.DataArray:
    da = rxr.open_rasterio(raster_path, masked=True).squeeze()  # (y,x)
    aoi = gpd.read_file(aoi_geojson)
    if da.rio.crs is None:
        raise ValueError(f"Raster has no CRS: {raster_path}")
    if aoi.crs is None:
        raise ValueError(f"AOI has no CRS: {aoi_geojson}")
    if str(aoi.crs) != str(da.rio.crs):
        aoi = aoi.to_crs(da.rio.crs)

    rb = da.rio.bounds()
    raster_poly = box(*rb)
    aoi_union = aoi.geometry.union_all()
    inter = raster_poly.intersection(aoi_union)
    if inter.is_empty:
        raise ValueError(f"No overlap between raster and AOI. Raster bounds={rb}, AOI bounds={tuple(aoi.total_bounds)}")

    # First try a normal clip dropping outside pixels
    try:
        clipped = da.rio.clip(aoi.geometry, aoi.crs, drop=True, invert=False)
    except NoDataInBounds:
        clipped = None

    # If result is empty or failed, crop to intersection bounds then mask w/ AOI without dropping grid
    if clipped is None or any(sz == 0 for sz in clipped.sizes.values()):
        minx, miny, maxx, maxy = inter.bounds
        # small padding to avoid zero-width due to rounding
        eps = max(abs(maxx-minx), abs(maxy-miny)) * 1e-9 or 1e-6
        cb = da.rio.clip_box(minx=minx-eps, miny=miny-eps, maxx=maxx+eps, maxy=maxy+eps)
        # Mask outside AOI but keep transform/grid
        try:
            cb_masked = cb.rio.clip(aoi.geometry, aoi.crs, drop=False, invert=False)
            clipped = cb_masked
        except NoDataInBounds:
            clipped = cb  # last resort: keep cropped box

    return clipped

def slope_from_dem(dem_path: Path) -> xr.DataArray:
    dem = rxr.open_rasterio(dem_path, masked=True).squeeze()
    transform = dem.rio.transform()
    px = transform.a; py = -transform.e
    dx = abs(px) if px else 1.0; dy = abs(py) if py else 1.0
    z = dem.values.astype("float64")
    if np.all(~np.isfinite(z)):
        return xr.full_like(dem, np.nan)
    dz_dy, dz_dx = np.gradient(z, dy, dx)
    slope_deg = np.degrees(np.arctan(np.sqrt(dz_dx**2 + dz_dy**2)))
    slope = xr.DataArray(slope_deg,
                         coords={dem.dims[0]: dem.coords[dem.dims[0]], dem.dims[1]: dem.coords[dem.dims[1]]},
                         dims=dem.dims, attrs=dem.attrs, name="Slope")
    if getattr(slope.rio, "crs", None) is None:
        slope = slope.rio.write_crs(dem.rio.crs)
    try:
        _ = slope.rio.transform(recalc=True)
    except Exception:
        slope = slope.rio.write_transform(dem.rio.transform(recalc=True))
    return slope

def reproject_match(src_da: xr.DataArray, match_da: xr.DataArray, resampling=Resampling.bilinear) -> xr.DataArray:
    out = src_da.rio.reproject_match(match_da, resampling=resampling)
    return _ensure_like(match_da, out, name=src_da.name or "reprojected")

def infiltration_score(awc_da: xr.DataArray, slope_da: xr.DataArray, lc_da: xr.DataArray | None,
                       w_awc: float, w_slp: float, w_lc: float) -> xr.DataArray:
    ref = awc_da
    awc_n = _normalize01(awc_da, method="percentile")
    slope_clamped = slope_da.clip(0, 30)
    slope_n = 1.0 - (slope_clamped / 30.0)
    slope_n = _ensure_like(ref, slope_n, name="SlopeIndex")
    if lc_da is None:
        lc_n = xr.full_like(ref, 0.5)
    else:
        lc_n = _normalize01(reproject_match(lc_da, ref, resampling=Resampling.nearest))

    w = np.array([w_awc, w_slp, w_lc], dtype="float64")
    if not np.isfinite(w).all() or w.sum() <= 0:
        w = np.array([1,1,1], dtype="float64")
    w = w / w.sum()

    comp = (w[0]*_ensure_like(ref, awc_n) +
            w[1]*_ensure_like(ref, slope_n) +
            w[2]*_ensure_like(ref, lc_n))
    score = (comp.clip(0,1) * 100.0).astype("float32")
    out = _ensure_like(ref, score, name="InfiltrationScore")
    return out
