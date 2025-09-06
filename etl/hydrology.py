from pathlib import Path
import numpy as np
import xarray as xr
import rioxarray as rxr
import geopandas as gpd
from shapely.geometry import box
from rioxarray.exceptions import NoDataInBounds
from rasterio.enums import Resampling

def clip_to_aoi(raster_path: Path, aoi_geojson: Path) -> xr.DataArray:
    da = rxr.open_rasterio(raster_path, masked=True).squeeze()
    aoi = gpd.read_file(aoi_geojson)
    if aoi.crs and str(aoi.crs) != str(da.rio.crs):
        aoi = aoi.to_crs(da.rio.crs)

    # Overlap guard
    rb = da.rio.bounds()
    inter = box(*rb).intersection(aoi.geometry.union_all())
    if inter.is_empty:
        raise ValueError(f"No overlap between raster and AOI. Raster bounds={rb}, AOI bounds={tuple(aoi.total_bounds)}")

    # Primary clip: keep grid, just mask
    try:
        clipped = da.rio.clip(aoi.geometry, aoi.crs, drop=False, invert=False)
    except NoDataInBounds:
        clipped = None

    # Fallback: crop to intersection bbox, pad ~0.5 pixel, then mask (still keep grid)
    if (clipped is None) or any(sz == 0 for sz in clipped.sizes.values()):
        minx, miny, maxx, maxy = inter.bounds
        rx, ry = da.rio.resolution()
        epsx, epsy = abs(rx) * 0.5, abs(ry) * 0.5
        cb = da.rio.clip_box(minx=minx - epsx, miny=miny - epsy, maxx=maxx + epsx, maxy=maxy + epsy)
        try:
            clipped = cb.rio.clip(aoi.geometry, aoi.crs, drop=False, invert=False)
        except NoDataInBounds:
            clipped = cb  # last resort: cropped window only

    return clipped

def reproject_match(da: xr.DataArray, ref: xr.DataArray, resampling=Resampling.nearest) -> xr.DataArray:
    return da.rio.reproject_match(ref, resampling=resampling)

def _to_da_like(ref: xr.DataArray, arr, name: str) -> xr.DataArray:
    if isinstance(arr, xr.DataArray) and tuple(arr.dims)==tuple(ref.dims):
        da = arr.copy()
    else:
        da = xr.DataArray(np.asarray(arr), coords={d: ref.coords[d] for d in ref.dims}, dims=ref.dims)
    da.name = name
    if getattr(da.rio, "crs", None) is None and getattr(ref.rio, "crs", None) is not None:
        da = da.rio.write_crs(ref.rio.crs)
    try:
        _ = da.rio.transform(recalc=True)
    except Exception:
        da = da.rio.write_transform(ref.rio.transform(recalc=True))
    return da

def slope_from_dem(dem_path: Path) -> xr.DataArray:
    da = rxr.open_rasterio(dem_path, masked=True).squeeze()
    try:
        slp = da.rio.slope()
    except Exception:
        y, x = np.gradient(da.values.astype("float32"))
        slp = np.degrees(np.arctan(np.sqrt(x*x + y*y)))
        slp = _to_da_like(da, slp, name="slope")
    slp = slp.clip(min=0)
    slp_norm = (slp.max() - slp) / (slp.max() - slp.min() + 1e-6)  # lower slope -> higher infiltration
    return _to_da_like(da, slp_norm, name="slope_norm")

def infiltration_score(awc_da: xr.DataArray, slope_da: xr.DataArray | None, lc_da: xr.DataArray | None,
                       w_awc: float, w_slope: float, w_lc: float) -> xr.DataArray:
    ref = awc_da
    A = awc_da.astype("float32")
    S = slope_da if slope_da is not None else xr.full_like(ref, 0.5, dtype="float32")
    if lc_da is None:
        L = xr.full_like(ref, 0.5, dtype="float32")
    else:
        L = lc_da

    # align to reference grid
    if getattr(S, "rio", None):
        try:
            S = S.rio.reproject_match(ref, resampling=Resampling.bilinear)
        except Exception:
            S = xr.full_like(ref, 0.5, dtype="float32")
    if getattr(L, "rio", None):
        try:
            L = L.rio.reproject_match(ref, resampling=Resampling.nearest)
        except Exception:
            L = xr.full_like(ref, 0.5, dtype="float32")

    score = (w_awc * A + w_slope * S + w_lc * L) / max(w_awc + w_slope + w_lc, 1e-6)
    score = score.clip(min=0.0, max=1.0)
    return _to_da_like(ref, score, name="InfiltrationScore")
