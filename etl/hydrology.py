"""
Hydrology ETL: clipping, reprojection, slope, and infiltration score.

Dependencies:
- numpy, xarray, rioxarray, rasterio
- shapely, geopandas
- etl.monitor.stage for timing/memory logs
"""

from __future__ import annotations
from pathlib import Path
from typing import Optional, Dict, Tuple

import numpy as np
import xarray as xr
import rioxarray  # noqa: F401  (registers .rio accessor)
import rasterio
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT

import geopandas as gpd
from shapely.geometry import box

try:
    # Shapely >=2
    from shapely import union_all as _union_all
except Exception:
    # Shapely <2
    from shapely.ops import unary_union as _union_all  # type: ignore

from etl.monitor import stage


# ---------------------------
# Paths / constants
# ---------------------------

def repo_root_from_file(__file__: str) -> Path:
    # etl/ -> project root
    return Path(__file__).resolve().parents[1]


_REPO_ROOT = repo_root_from_file(__file__)
DATA_DIR = _REPO_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROC_DIR = DATA_DIR / "processed"
HYDRO_DIR = PROC_DIR / "hydrology"
AOI_DIR = PROC_DIR / "aoi"
HYDRO_DIR.mkdir(parents=True, exist_ok=True)
AOI_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------
# Utilities
# ---------------------------

def _ensure_da(da_or_ds: xr.DataArray | xr.Dataset, band: Optional[str] = None) -> xr.DataArray:
    if isinstance(da_or_ds, xr.Dataset):
        if band is None:
            # pick the first data variable
            band = list(da_or_ds.data_vars)[0]
        return da_or_ds[band]
    return da_or_ds


def open_raster(path: str | Path) -> xr.DataArray:
    """Open a raster as a DataArray with a CRS via rioxarray."""
    path = Path(path)
    da = xr.open_dataarray(path)
    # In case it opened as Dataset
    da = _ensure_da(da)
    da = da.rio.write_crs(da.rio.crs or rasterio.open(path).crs, inplace=True)
    return da


def save_raster(da: xr.DataArray, path: str | Path, dtype: Optional[str] = None, nodata: Optional[float] = None) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    if dtype:
        da = da.astype(dtype)
    if nodata is not None:
        # rioxarray write_crs/write_nodata for explicit NODATA
        da = da.rio.write_nodata(nodata, encoded=True)
    da.rio.to_raster(p)
    return p


# ---------------------------
# Spatial transforms
# ---------------------------

def reproject_match(source: xr.DataArray, match: xr.DataArray, resampling: Resampling = Resampling.bilinear) -> xr.DataArray:
    """Reproject `source` to `match` grid (crs, transform, shape)."""
    src = source
    src = src.rio.reproject_match(match, resampling=resampling)
    return src


def clip_to_aoi(da: xr.DataArray, aoi: gpd.GeoDataFrame) -> xr.DataArray:
    """
    Clip a raster to an AOI GeoDataFrame. Handles CRS differences.
    """
    if getattr(aoi, "crs", None) is None:
        raise ValueError("AOI has no CRS; expected a projected CRS like EPSG:25830")

    if da.rio.crs is None:
        raise ValueError("Raster has no CRS.")

    # Reproject AOI to raster CRS if needed
    if aoi.crs != da.rio.crs:
        aoi_r = aoi.to_crs(da.rio.crs)
    else:
        aoi_r = aoi

    # Check overlap
    rb = da.rio.bounds()
    inter = box(*rb).intersection(_union_all(list(aoi_r.geometry)))
    if inter.is_empty:
        raise ValueError(
            f"No overlap between raster and AOI. Raster bounds={rb}, AOI bounds={tuple(aoi_r.total_bounds)}"
        )

    clipped = da.rio.clip(aoi_r.geometry, aoi_r.crs, drop=True)
    clipped = clipped.where(np.isfinite(clipped), np.nan)
    return clipped


def slope_from_dem(dem_da: xr.DataArray, pixel_size: Optional[Tuple[float, float]] = None) -> xr.DataArray:
    """
    Compute slope (in degrees) from DEM.
    Uses simple 3x3 Sobel-like kernel in projected units.
    """
    dem = dem_da
    # Estimate pixel size in meters from transform if not given
    if pixel_size is None:
        try:
            t = dem.rio.transform()
            dx = abs(t.a)
            dy = abs(t.e)
        except Exception:
            dx = dy = 1.0
    else:
        dx, dy = pixel_size

    # Pad edges with nan
    z = dem.values
    if z.ndim == 3:
        # (band, y, x) -> take first band
        z = z[0]
    z = z.astype("float64")

    # Build gradients (finite differences)
    gy, gx = np.gradient(z, dy, dx)
    slope_rad = np.arctan(np.hypot(gx, gy))
    slope_deg = np.degrees(slope_rad)

    out = xr.DataArray(
        slope_deg,
        coords=dem.coords,
        dims=dem.dims,
        attrs=dem.attrs
    )
    out.rio.write_crs(dem.rio.crs, inplace=True)
    return out


def infiltration_score(awc_da: xr.DataArray, slope_da: xr.DataArray, lc_factor: Optional[xr.DataArray] = None) -> xr.DataArray:
    """
    Combine AWC (available water capacity), slope, and optional land-cover factor into an "infiltration score".

    Example scaling (simple, for MVP):
      score = sigmoid( z_awc - z_slope + lc_term )

    where:
      z_awc, z_slope are standardized (z-scores) within the AOI,
      lc_term is 0 if not provided, otherwise scaled to [-0.5, +0.5].
    """
    # Match grids
    slope_da = reproject_match(slope_da, awc_da, resampling=Resampling.bilinear)
    if lc_factor is not None:
        lc_factor = reproject_match(lc_factor, awc_da, resampling=Resampling.nearest)

    # Convert to numpy, ignore nans for stats
    def zscore(arr):
        m = np.nanmean(arr)
        s = np.nanstd(arr)
        if not np.isfinite(s) or s == 0:
            return np.zeros_like(arr)
        return (arr - m) / s

    awc = awc_da.values.astype("float64")
    slope = slope_da.values.astype("float64")

    z_awc = zscore(awc)
    z_slope = zscore(slope)

    lc_term = 0.0
    if lc_factor is not None:
        lc = lc_factor.values.astype("float64")
        # Scale to [-0.5, 0.5] based on percentiles
        lo, hi = np.nanpercentile(lc, [5, 95])
        rng = (hi - lo) if hi > lo else 1.0
        lc_term = ((lc - lo) / rng) - 0.5

    z = z_awc - z_slope + lc_term

    # Sigmoid to [0,1]
    score = 1.0 / (1.0 + np.exp(-z))
    out = xr.DataArray(
        score,
        coords=awc_da.coords,
        dims=awc_da.dims,
        attrs=awc_da.attrs
    )
    out.rio.write_crs(awc_da.rio.crs, inplace=True)
    return out


# ---------------------------
# Pipeline orchestrator (optional)
# ---------------------------

def run_pipeline(
    awc_path: str | Path,
    dem_path: str | Path,
    aoi_gpkg: str | Path,
    out_dir: str | Path | None = None,
    lc_path: Optional[str | Path] = None,
) -> Dict[str, Path]:
    """
    Minimal orchestrator that:
      - opens AWC & DEM
      - clips to AOI
      - computes slope
      - computes infiltration
      - saves outputs under data/processed/hydrology/
    Returns dict of output paths.
    """
    out_dir = Path(out_dir or HYDRO_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)

    awc = open_raster(awc_path)
    dem = open_raster(dem_path)
    aoi = gpd.read_file(aoi_gpkg)
    with stage("clip AWC"):
        awc_c = clip_to_aoi(awc, aoi)
    with stage("clip DEM"):
        dem_c = clip_to_aoi(dem, aoi)
    with stage("slope"):
        slope = slope_from_dem(dem_c)
    if lc_path:
        lc = open_raster(lc_path)
        with stage("clip LC"):
            lc_c = clip_to_aoi(lc, aoi)
    else:
        lc_c = None
    with stage("infiltration"):
        infil = infiltration_score(awc_c, slope, lc_c)

    # Save
    awc_out = save_raster(awc_c, out_dir / "AWC_clip.tif", dtype="float32", nodata=np.nan)
    dem_out = save_raster(dem_c, out_dir / "DEM_clip.tif", dtype="float32", nodata=np.nan)
    slope_out = save_raster(slope, out_dir / "Slope_deg.tif", dtype="float32", nodata=np.nan)
    infil_out = save_raster(infil, out_dir / "InfiltrationScore.tif", dtype="float32", nodata=np.nan)

    return {
        "awc": awc_out,
        "dem": dem_out,
        "slope": slope_out,
        "infiltration": infil_out,
    }


# ---------------------------
# Sampling helper (for map inspector)
# ---------------------------

def sample_rasters_at(lon: float, lat: float, layers: Dict[str, str]) -> Dict[str, float]:
    """
    Sample multiple rasters by lon/lat (EPSG:4326).
    layers: {"awc": "/path.tif", "slope": "...", "infiltration": "..."}
    """
    out: Dict[str, float] = {}
    for name, p in layers.items():
        path = Path(p)
        if not path.exists():
            out[name] = float("nan")
            continue
        with rasterio.open(path) as src:
            with WarpedVRT(src, crs="EPSG:4326", resampling=Resampling.bilinear) as vrt:
                v = list(vrt.sample([(lon, lat)]))[0][0]
        out[name] = float(v) if v is not None and np.isfinite(v) else float("nan")
    return out
