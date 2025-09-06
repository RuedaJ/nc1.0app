from pathlib import Path
import rasterio as rio
import numpy as np
import rioxarray as rxr
import xarray as xr
from rasterio.enums import Resampling
import geopandas as gpd

def reproject_match(src_da: xr.DataArray, ref_da: xr.DataArray, resampling=Resampling.bilinear) -> xr.DataArray:
    return src_da.rio.reproject_match(ref_da, resampling=resampling)

def slope_from_dem(dem_path: Path) -> xr.DataArray:
    da = rxr.open_rasterio(dem_path, masked=True).squeeze()
    x, y = np.gradient(da.values, da.rio.resolution()[0], da.rio.resolution()[1])
    slope_rad = np.arctan(np.hypot(x, y))
    slope_deg = np.degrees(slope_rad)
    out = xr.DataArray(slope_deg, coords=da.coords, dims=da.dims, attrs=da.attrs, name="slope")
    out.rio.write_nodata(np.nan, inplace=True)
    return out

def minmax_norm(da: xr.DataArray, clip=(0.0, 99.0)) -> xr.DataArray:
    lo, hi = np.nanpercentile(da.values, clip)
    scaled = (da - lo) / max(1e-9, (hi - lo))
    return scaled.clip(0,1)

def infiltration_score(awc_da: xr.DataArray, slope_da: xr.DataArray, lc_da: xr.DataArray|None,
                       w_awc: float, w_slope: float, w_lc: float) -> xr.DataArray:
    awc_n = minmax_norm(awc_da)
    slope_n = minmax_norm(slope_da)
    slope_inv = 1 - slope_n
    if lc_da is not None:
        lc_n = minmax_norm(lc_da)
    else:
        lc_n = xr.full_like(awc_n, fill_value=0.5)

    score = w_awc*awc_n + w_slope*slope_inv + w_lc*lc_n
    score = (score*100).clip(0,100)
    score.name = "InfiltrationScore"
    score.rio.write_nodata(np.nan, inplace=True)
    return score

def clip_to_aoi(raster_path: Path, aoi_geojson: Path) -> xr.DataArray:
    da = rxr.open_rasterio(raster_path, masked=True).squeeze()
    aoi = gpd.read_file(aoi_geojson)
    if aoi.crs and str(aoi.crs) != str(da.rio.crs):
        aoi = aoi.to_crs(da.rio.crs)
    clipped = da.rio.clip(aoi.geometry, aoi.crs, drop=True, invert=False)
    return clipped
