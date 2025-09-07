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
    if da.ndim == 3:
        da = da.squeeze(drop=True)
    return da


def clip_to_aoi(raster_path: str | Path, aoi_geojson: str | Path) -> xr.DataArray:
    da = _ensure_2d(rxr.open_rasterio(raster_path, masked=True))
    aoi = gpd.read_file(aoi_geojson)
    if aoi.empty:
        raise ValueError("AOI is empty.")
    if aoi.crs is not None and da.rio.crs is not None and str(aoi.crs) != str(da.rio.crs):
        aoi_r = aoi.to_crs(da.rio.crs)
    else:
        aoi_r = aoi

    rb = da.rio.bounds()
    inter = box(*rb).intersection(_union_all(list(aoi_r.geometry)))
    if inter.is_empty:
        raise ValueError(
            f"No overlap between raster and AOI. "
            f"Raster bounds={rb}, AOI bounds={tuple(aoi_r.total_bounds)}"
        )

    clipped = da.rio.clip(aoi_r.geometry, aoi_r.crs, drop=True)
    clipped = clipped.where(np.isfinite(clipped), np.nan)
    if getattr(clipped.rio, "crs", None) is None and getattr(da.rio, "crs", None) is not None:
        clipped = clipped.rio.write_crs(da.rio.crs)
    try:
        _ = clipped.rio.transform(recalc=True)
    except Exception:
        clipped = clipped.rio.write_transform(da.rio.transform(recalc=True))
    return clipped


def slope_from_dem(dem_path: str | Path) -> xr.DataArray:
    dem = _ensure_2d(rxr.open_rasterio(dem_path, masked=True)).where(lambda x: np.isfinite(x), np.nan)
    resx, resy = dem.rio.resolution()
    arr = dem.values.astype("float32")
    with np.errstate(invalid="ignore", divide="ignore"):
        gy, gx = np.gradient(arr, abs(resy), abs(resx))
        slope_rad = np.arctan(np.hypot(gx, gy))
        slope_deg = np.degrees(slope_rad)
    slope = xr.DataArray(slope_deg, coords=dem.coords, dims=dem.dims, attrs=dem.attrs, name="slope_deg")
    if dem.rio.crs is not None:
        slope = slope.rio.write_crs(dem.rio.crs, inplace=False)
        slope = slope.rio.write_transform(dem.rio.transform(recalc=True), inplace=False)
    return slope


def reproject_match(src: xr.DataArray, template: xr.DataArray,
                    resampling: Resampling = Resampling.nearest) -> xr.DataArray:
    src = _ensure_2d(src); template = _ensure_2d(template)
    if src.rio.crs is None:
        raise ValueError("Source raster has no CRS; cannot reproject_match.")
    if template.rio.crs is None:
        raise ValueError("Template raster has no CRS; cannot reproject_match.")
    out = src.rio.reproject_match(template, resampling=resampling)
    return out.where(np.isfinite(out), np.nan)


def _norm_0_1(da: xr.DataArray, invert: bool = False) -> xr.DataArray:
    arr = da.values.astype("float32")
    finite = np.isfinite(arr)
    if finite.any():
        vmin = float(np.nanmin(arr[finite])); vmax = float(np.nanmax(arr[finite]))
        norm = (arr - vmin) / (vmax - vmin) if vmax > vmin else np.zeros_like(arr, dtype="float32")
    else:
        norm = np.zeros_like(arr, dtype="float32")
    if invert:
        norm = 1.0 - norm
    out = xr.DataArray(norm, coords=da.coords, dims=da.dims, attrs=da.attrs)
    if da.rio.crs is not None:
        out = out.rio.write_crs(da.rio.crs, inplace=False)
        out = out.rio.write_transform(da.rio.transform(recalc=True), inplace=False)
    return out


def _same_grid(a: xr.DataArray, b: xr.DataArray) -> bool:
    try:
        same_shape = a.sizes == b.sizes
        same_crs = str(a.rio.crs) == str(b.rio.crs)
        same_tx = tuple(a.rio.transform()) == tuple(b.rio.transform())
        return bool(same_shape and same_crs and same_tx)
    except Exception:
        return False


def infiltration_score(awc_da: xr.DataArray,
                       slope_da: xr.DataArray,
                       lc_da: xr.DataArray | None,
                       w_awc: float, w_slope: float, w_lc: float) -> xr.DataArray:
    """
    Compute infiltration score in [0,1]:
    - Higher AWC -> higher score
    - Flatter slope -> higher score (slope inverted)
    - Land cover optional (defaults to 0.5)
    All inputs are aligned to the AWC grid before combining.
    """
    # Use AWC grid as the template
    tmpl = _ensure_2d(awc_da)

    # Align slope to AWC grid
    slope_aligned = _ensure_2d(slope_da)
    if not _same_grid(slope_aligned, tmpl):
        slope_aligned = reproject_match(slope_aligned, tmpl, resampling=Resampling.bilinear)

    # Align LC (if any) to AWC grid
    if lc_da is None:
        lc_aligned = xr.full_like(tmpl, 0.5)
    else:
        lc_aligned = _ensure_2d(lc_da)
        if not _same_grid(lc_aligned, tmpl):
            lc_aligned = reproject_match(lc_aligned, tmpl, resampling=Resampling.nearest)

    # Normalize on the aligned grids
    awc_n   = _norm_0_1(tmpl, invert=False)
    slope_n = _norm_0_1(slope_aligned, invert=True)  # flatter is better
    lc_n    = _norm_0_1(lc_aligned, invert=False) if lc_da is not None else lc_aligned

    weights = np.array([w_awc, w_slope, w_lc], dtype="float32")
    sw = float(weights.sum()) if float(weights.sum()) > 0 else 1.0
    score_vals = (w_awc * awc_n.values + w_slope * slope_n.values + w_lc * lc_n.values) / sw
    score_vals = np.clip(score_vals, 0.0, 1.0).astype("float32")

    score = xr.DataArray(score_vals, coords=tmpl.coords, dims=tmpl.dims, attrs=tmpl.attrs, name="InfiltrationScore")
    if tmpl.rio.crs is not None:
        score = score.rio.write_crs(tmpl.rio.crs, inplace=False)
        score = score.rio.write_transform(tmpl.rio.transform(recalc=True), inplace=False)
    return score
