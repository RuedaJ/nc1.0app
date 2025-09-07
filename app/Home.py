import sys, os
_REPO_ROOT = os.path.dirname(os.path.dirname(__file__))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import streamlit as st
import yaml
from pathlib import Path
import numpy as np
import xarray as xr

from state import get_state, reset
from etl.validators import valid_lat, valid_lon, looks_like_drive_url, validate_geotiff, sniff_signature
from etl.drive import parse_drive_url, download_drive_any, download_to
from etl.geometry import aoi_from_latlon, save_geojson
from etl.hydrology import clip_to_aoi, slope_from_dem, reproject_match, infiltration_score
from rasterio.enums import Resampling
import rasterio as rio
import geopandas as gpd
from shapely.geometry import box

st.set_page_config(page_title="sustai-geo-app", layout="wide")

try:
    cap = st.runtime.get_option("server.maxUploadSize")
except Exception:
    cap = st.config.get_option("server.maxUploadSize")
st.caption(f"Max upload size (MB) in effect: {cap}")


cfg = yaml.safe_load(Path("configs/default.yaml").read_text())
S = get_state()

st.sidebar.markdown("### 🌍 sustai-geo-app")
lang = st.sidebar.selectbox("Language/Idioma", ["EN","ES"],
    index=0 if cfg["app"]["language_default"]=="EN" else 1, key="lang_selector")

st.sidebar.markdown("---")
st.sidebar.markdown("**Site Location**")
lat = st.sidebar.text_input("Latitude", key="lat_input", placeholder="e.g., 41.6488")
lon = st.sidebar.text_input("Longitude", key="lon_input", placeholder="e.g., -0.8891")
st.sidebar.button("📍 Use map click", key="map_click_btn")

st.sidebar.markdown("**Analysis Area**")
radius_km = st.sidebar.slider("Analysis Radius (km)", 1, cfg["app"]["max_radius_km"],
    S.get("radius_km", cfg["app"]["default_radius_km"]), key="radius_slider")

st.sidebar.markdown("**Upload rasters (preferred)**")
awc_up = st.sidebar.file_uploader("AWC GeoTIFF (*.tif, *.tiff, *.vrt)", type=["tif","tiff","vrt"], key="awc_upl")
dem_up = st.sidebar.file_uploader("DEM GeoTIFF (*.tif, *.tiff, *.vrt)", type=["tif","tiff","vrt"], key="dem_upl")

st.sidebar.markdown("**…or provide links/paths (fallback)**")
st.sidebar.text_input("AWC Raster link or path", value=S.get("awc_link") or "", key="awc_link_input")
st.sidebar.text_input("DEM (PNOA) link or path", value=S.get("dem_link") or "", key="dem_link_input")
st.sidebar.text_input("Land Cover (optional link or path)", value=S.get("lc_link") or "", key="lc_link_input")

with st.sidebar.expander("Advanced Settings"):
    st.slider("AWC Weight", 0.0, 1.0, cfg["weights"]["water"]["awc"], 0.05, key="awc_weight")
    st.slider("Slope Weight", 0.0, 1.0, cfg["weights"]["water"]["slope"], 0.05, key="slope_weight")
    st.slider("Land Cover Weight", 0.0, 1.0, cfg["weights"]["water"]["landcover"], 0.05, key="lc_weight")

run = st.sidebar.button("🚀 Run Analysis", use_container_width=True, key="run_analysis")
if st.sidebar.button("🗑️ Clear Cache", key="clear_cache"):
    reset(); st.sidebar.success("Cache cleared and session reset.")

def _resolve_raster_file(pathlike: Path, preferred_names=(
    "awc.tif","AWC.tif","awc.vrt",
    "dem.tif","DEM.tif","dem.vrt"
)) -> Path:
    p = Path(pathlike)
    if p.is_file():
        return p
    if p.is_dir():
        for name in preferred_names:
            cand = p / name
            if cand.is_file(): return cand
        for ext in (".tif",".tiff",".vrt",".img"):
            cands = sorted(p.glob(f"*{ext}"))
            if cands: return cands[0]
    return p

def _save_upload(upl, out: Path):
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "wb") as f:
        f.write(upl.read())
    return out

def _reject_part(p: Path, label: str):
    if p.suffix.lower() == ".part" or p.name.endswith(".part"):
        st.error(f"{label} looks incomplete ({p.name}). Click '🗑️ Clear Cache' and upload again.")
        st.stop()

def _ensure_like(ref: xr.DataArray, out, name="InfiltrationScore"):
    if isinstance(out, xr.DataArray) and tuple(out.dims)==tuple(ref.dims) and all(d in out.coords for d in ref.dims):
        da = out.copy(); da.name = name
    else:
        da = xr.DataArray(np.asarray(out),
                          coords={ref.dims[0]: ref.coords[ref.dims[0]], ref.dims[1]: ref.coords[ref.dims[1]]},
                          dims=ref.dims, attrs=ref.attrs, name=name)
    if getattr(da.rio, "crs", None) is None and getattr(ref.rio, "crs", None) is not None:
        da = da.rio.write_crs(ref.rio.crs)
    try:
        _ = da.rio.transform(recalc=True)
    except Exception:
        da = da.rio.write_transform(ref.rio.transform(recalc=True))
    return da

st.title("Welcome 👋")
st.write("Upload local rasters (recommended). You can still use links/paths as fallback.")
st.info("After a successful run, switch to the Mapa page to see the Infiltration layer.")

if run:
    errors = []
    if not valid_lat(lat) or not valid_lon(lon):
        errors.append("Invalid coordinates. Please enter valid latitude and longitude.")
    _awc_link = st.session_state.get("awc_link_input") or ""
    _dem_link = st.session_state.get("dem_link_input") or ""
    _lc_link  = st.session_state.get("lc_link_input") or ""

    if errors:
        for e in errors: st.error(e); S["errors"]=errors
        st.stop()

    S["coords"]=(float(lat), float(lon)); S["radius_km"]=int(radius_km)
    S["awc_link"]=_awc_link; S["dem_link"]=_dem_link; S["lc_link"]=(_lc_link or None)

    # AOI
    aoi_gdf = aoi_from_latlon(float(lat), float(lon), float(radius_km), cfg["app"]["crs"])
    aoi_path = Path(cfg["paths"]["processed"]) / "site_aoi" / f"site_buffer_{radius_km}km.geojson"
    save_geojson(aoi_gdf, aoi_path)
    S["aoi_geojson"]=str(aoi_path)

    raw_dir = Path(cfg["paths"]["raw"])

    # AWC: prefer local upload, then link/path fallback
    if awc_up is not None:
        awc_raw = _save_upload(awc_up, raw_dir/"awc"/"awc.tif")
    else:
        if Path(_awc_link).exists():
            awc_raw = Path(_awc_link)
        elif parse_drive_url(_awc_link)[0] in ("file","folder"):
            awc_raw = download_drive_any(_awc_link, raw_dir/"awc")
        elif _awc_link:
            awc_raw = download_to(raw_dir/"awc"/"awc.tif", _awc_link)
        else:
            st.error("Please upload AWC or provide a link/path."); st.stop()
    awc_raw = _resolve_raster_file(awc_raw)
    _reject_part(awc_raw, "AWC")

    # DEM: prefer local upload, then link/path fallback
    if dem_up is not None:
        dem_raw = _save_upload(dem_up, raw_dir/"dem"/"dem.tif")
    else:
        if Path(_dem_link).exists():
            dem_raw = Path(_dem_link)
        elif parse_drive_url(_dem_link)[0] in ("file","folder"):
            dem_raw = download_drive_any(_dem_link, raw_dir/"dem")
        elif _dem_link:
            dem_raw = download_to(raw_dir/"dem"/"dem.tif", _dem_link)
        else:
            st.error("Please upload DEM or provide a link/path."); st.stop()
    dem_raw = _resolve_raster_file(dem_raw)
    _reject_part(dem_raw, "DEM")

    # Diagnostics: signatures
    awc_sig = sniff_signature(str(awc_raw)) if awc_raw.exists() and awc_raw.is_file() else ("DIR" if awc_raw.exists() and awc_raw.is_dir() else "MISSING")
    dem_sig = sniff_signature(str(dem_raw)) if dem_raw.exists() and dem_raw.is_file() else ("DIR" if dem_raw.exists() and dem_raw.is_dir() else "MISSING")
    st.caption(f"AWC path: {awc_raw} • size={(awc_raw.stat().st_size if awc_raw.exists() and awc_raw.is_file() else 'NA')} • sig={awc_sig}")
    st.caption(f"DEM path: {dem_raw} • size={(dem_raw.stat().st_size if dem_raw.exists() and dem_raw.is_file() else 'NA')} • sig={dem_sig}")

    # Validate files
    ok_awc, msg_awc = (False,"Path is a directory") if awc_raw.is_dir() else validate_geotiff(str(awc_raw))
    ok_dem, msg_dem = (False,"Path is a directory") if dem_raw.is_dir() else validate_geotiff(str(dem_raw))
    if not ok_awc: st.error(f"AWC invalid: {msg_awc}"); st.stop()
    if not ok_dem: st.error(f"DEM invalid: {msg_dem}"); st.stop()

    S["awc_path"], S["dem_path"]=str(awc_raw), str(dem_raw)
    S["awc_loaded"], S["dem_loaded"]=True, True

    # Overlap diagnostics for DEM & AWC
    aoi = gpd.read_file(aoi_path)

    with rio.open(str(dem_raw)) as ds:
        dem_bounds = ds.bounds; dem_crs = ds.crs
    aoi_r_dem = aoi.to_crs(dem_crs) if str(aoi.crs) != str(dem_crs) else aoi
    dem_overlap = box(*dem_bounds).intersects(aoi_r_dem.geometry.union_all())
    st.caption(f"DEM bounds: {tuple(round(v,5) for v in dem_bounds)} | AOI overlaps DEM: {dem_overlap}")

    with rio.open(str(awc_raw)) as ds2:
        awc_bounds = ds2.bounds; awc_crs = ds2.crs
    aoi_r_awc = aoi.to_crs(awc_crs) if str(aoi.crs) != str(awc_crs) else aoi
    awc_overlap = box(*awc_bounds).intersects(aoi_r_awc.geometry.union_all())
    st.caption(f"AWC bounds: {tuple(round(v,5) for v in awc_bounds)} | AOI overlaps AWC: {awc_overlap}")
    if not awc_overlap:
        st.error("AOI does not overlap AWC. Adjust coordinates/radius or select an AWC covering the AOI.")
        st.stop()

    # Clip & compute
    awc_da_clip = clip_to_aoi(awc_raw, aoi_path)
    st.caption(f"AWC clip sizes: {[awc_da_clip.sizes[d] for d in awc_da_clip.dims]}")
    if any(sz == 0 for sz in awc_da_clip.sizes.values()):
        st.error("AWC clip resulted in empty array. Check that AWC overlaps the AOI or use a larger radius.")
        st.stop()

    try:
        dem_da_clip = clip_to_aoi(dem_raw, aoi_path)
        if any(sz == 0 for sz in dem_da_clip.sizes.values()):
            raise ValueError("DEM clip returned empty grid")
        slope_da    = slope_from_dem(dem_raw)
        slope_da    = reproject_match(slope_da, dem_da_clip, resampling=Resampling.bilinear)
        dem_ok = True
    except Exception as e:
        st.warning(f"DEM could not be clipped to AOI ({e}). Falling back to neutral slope (0.5).")
        slope_da = xr.full_like(awc_da_clip, 0.5)
        dem_ok = False
        S["dem_loaded"] = False

    # Compute score
    w_awc=float(st.session_state["awc_weight"]); w_slp=float(st.session_state["slope_weight"]); w_lc=float(st.session_state.get("lc_weight", cfg["weights"]["water"]["landcover"]))
    score = infiltration_score(awc_da_clip, slope_da, None, w_awc, w_slp, w_lc)
    score = _ensure_like(awc_da_clip, score, name="InfiltrationScore")

    # --- ensure CRS & transform before save ---
    score = score.where(np.isfinite(score), np.nan)
    score = score.rio.write_nodata(np.nan, inplace=False)
    ref_crs = awc_da_clip.rio.crs
    ref_transform = awc_da_clip.rio.transform(recalc=True)
    score = score.rio.write_crs(ref_crs, inplace=False)
    score = score.rio.write_transform(ref_transform, inplace=False)

    # Debug captions
    y_dim, x_dim = score.dims
    st.caption(f"Score dims={score.dims}, sizes=({score.sizes[y_dim]}, {score.sizes[x_dim]})")
    st.caption(f"Score CRS set? {score.rio.crs is not None}")
    try:
        _ = score.rio.transform(recalc=False)
        st.caption("Score has cached transform ✅")
    except Exception:
        st.caption("Score transform missing; will use explicit write/fallback.")

    if score.sizes[y_dim] == 0 or score.sizes[x_dim] == 0:
        st.error("InfiltrationScore grid is empty after processing.")
        st.stop()

    # Save with rioxarray; fallback to rasterio
    out_dir = Path(cfg["paths"]["processed"]) / "hydrology"; out_dir.mkdir(parents=True, exist_ok=True)
    out_tif = out_dir / "InfiltrationScore.tif"
    try:
        score.rio.to_raster(out_tif, compress="LZW", recalc_transform=False)
    except Exception as e:
        h, w = int(score.sizes[y_dim]), int(score.sizes[x_dim])
        profile = {
            "driver": "GTiff",
            "height": h,
            "width":  w,
            "count":  1,
            "dtype":  "float32",
            "crs":    ref_crs,
            "transform": ref_transform,
            "compress": "LZW",
            "nodata": np.nan,
        }
        with rio.open(out_tif, "w", **profile) as dst:
            dst.write(score.values.astype("float32"), 1)
        st.warning(f"Saved via rasterio fallback due to rioxarray error: {e}")

    S["infiltration_path"]=str(out_tif)
    finite_ratio = float(np.isfinite(score.values).mean())
    st.caption(f"InfiltrationScore finite ratio: {finite_ratio:.3f}")
    st.success("✅ Analysis completed. Open **Mapa** to view the Infiltration layer." if dem_ok else "✅ Analysis completed (neutral slope). Open **Mapa**.")
