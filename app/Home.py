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

from etl.validators import (
    valid_lat, valid_lon, looks_like_drive_url,
    validate_geotiff, sniff_signature,
)
from etl.drive import extract_id, download_drive_file, download_to
from etl.geometry import aoi_from_latlon, save_geojson
from etl.hydrology import (
    clip_to_aoi, slope_from_dem, reproject_match, infiltration_score
)
from rasterio.enums import Resampling
import rasterio as rio
import geopandas as gpd
from shapely.geometry import box

st.set_page_config(page_title="sustai-geo-app", layout="wide")

cfg = yaml.safe_load(Path("configs/default.yaml").read_text())
S = get_state()

st.sidebar.markdown("### 🌍 sustai-geo-app")
lang = st.sidebar.selectbox(
    "Language/Idioma", ["EN","ES"],
    index=0 if cfg["app"]["language_default"]=="EN" else 1,
    key="lang_selector"
)

st.sidebar.markdown("---")
st.sidebar.markdown("**Site Location**")
lat = st.sidebar.text_input("Latitude", key="lat_input", placeholder="e.g., 41.6488")
lon = st.sidebar.text_input("Longitude", key="lon_input", placeholder="e.g., -0.8891")
st.sidebar.button("📍 Use map click", key="map_click_btn")

st.sidebar.markdown("**Analysis Area**")
radius_km = st.sidebar.slider("Analysis Radius (km)", 1, cfg["app"]["max_radius_km"], S.get("radius_km", cfg["app"]["default_radius_km"]), key="radius_slider")

st.sidebar.markdown("**Data Sources (Google Drive or local path)**")
st.sidebar.text_input("AWC Raster Link *", value=S.get("awc_link") or "", key="awc_link_input", help="Drive link OR local path to a .tif/.vrt (or a folder containing one)")
st.sidebar.text_input("DEM (PNOA) Link *", value=S.get("dem_link") or "", key="dem_link_input", help="Drive link OR local path to a .tif/.vrt (or a folder containing one)")
st.sidebar.text_input("Land Cover (optional)", value=S.get("lc_link") or "", key="lc_link_input")

with st.sidebar.expander("Advanced Settings"):
    st.slider("AWC Weight", 0.0, 1.0, cfg["weights"]["water"]["awc"], 0.05, key="awc_weight")
    st.slider("Slope Weight", 0.0, 1.0, cfg["weights"]["water"]["slope"], 0.05, key="slope_weight")
    st.slider("Land Cover Weight", 0.0, 1.0, cfg["weights"]["water"]["landcover"], 0.05, key="lc_weight")

run = st.sidebar.button("🚀 Run Analysis", use_container_width=True, key="run_analysis")
if st.sidebar.button("🗑️ Clear Cache", key="clear_cache"):
    reset(); st.sidebar.success("Cache cleared and session reset.")

status = []
status.append("AWC ✓" if S.get("awc_loaded") else "AWC ▫")
status.append("DEM ✓" if S.get("dem_loaded") else "DEM ▫")
status.append("APIs ✓" if S.get("api_ok") else "APIs ▫")
st.sidebar.markdown("---")
st.sidebar.caption("Status: " + "  ".join(status))

# ---- helpers to resolve directories -> raster files ----
RASTER_EXTS = (".tif", ".tiff", ".vrt", ".img")  # extend as needed

def _resolve_raster_file(pathlike: Path, preferred_names=("awc.tif","dem.tif")) -> Path:
    p = Path(pathlike)
    if p.is_file():
        return p
    if p.is_dir():
        # try preferred names
        for name in preferred_names:
            cand = p / name
            if cand.exists() and cand.is_file():
                return cand
        # else first georaster-looking file
        for ext in RASTER_EXTS:
            cands = sorted(p.glob(f"*{ext}"))
            if cands:
                return cands[0]
    return p  # could still be invalid; validator will catch

def _dir_hint(p: Path, max_items=5):
    try:
        if not p.is_dir():
            return None
        items = sorted([q.name for q in p.iterdir()])[:max_items]
        return ", ".join(items) if items else "(empty)"
    except Exception:
        return None

# ───────────────────── Main ─────────────────────
st.title("Welcome 👋")
st.write("Use the tabs to navigate: **Mapa**, **Zonificación**, **Dashboard**, **Informe**.")
st.info("Enter coordinates and Google Drive links (or local paths) in the sidebar, then click **Run Analysis**. Data will be cached locally.")

if run:
    errors = []
    if not valid_lat(lat) or not valid_lon(lon):
        errors.append("Invalid coordinates. Please enter valid latitude and longitude.")
    _awc = st.session_state.get("awc_link_input") or ""
    _dem = st.session_state.get("dem_link_input") or ""
    _lc  = st.session_state.get("lc_link_input") or ""

    if cfg["inputs"]["require_drive_links"]:
        from pathlib import Path as _P
        if not (_P(_awc).exists() or looks_like_drive_url(_awc)): errors.append("AWC: Invalid Google Drive link.")
        if not (_P(_dem).exists() or looks_like_drive_url(_dem)): errors.append("DEM: Invalid Google Drive link.")

    if errors:
        for e in errors: st.error(e); S["errors"]=errors
        st.stop()

    S["coords"]=(float(lat), float(lon)); S["radius_km"]=int(radius_km)
    S["awc_link"]=_awc; S["dem_link"]=_dem; S["lc_link"]=(_lc or None)

    # AOI
    aoi_gdf = aoi_from_latlon(float(lat), float(lon), float(radius_km), cfg["app"]["crs"])
    aoi_path = Path(cfg["paths"]["processed"]) / "site_aoi" / f"site_buffer_{radius_km}km.geojson"
    save_geojson(aoi_gdf, aoi_path)
    S["aoi_geojson"]=str(aoi_path)

    raw_dir = Path(cfg["paths"]["raw"])

    # AWC
    from etl.drive import download_drive_any, parse_drive_url
    if Path(_awc).exists():
        awc_raw = Path(_awc)
    elif parse_drive_url(_awc)[0] in ("file","folder"):
        awc_raw = download_drive_any(_awc, raw_dir/"awc")
    else:
        awc_raw = download_to(raw_dir/"awc"/"awc.tif", _awc)
    awc_raw = _resolve_raster_file(awc_raw, ("awc.tif","AWC.tif","awc.vrt"))

    # DEM
    if Path(_dem).exists():
        dem_raw = Path(_dem)
    elif parse_drive_url(_dem)[0] in ("file","folder"):
        dem_raw = download_drive_any(_dem, raw_dir/"dem")
    else:
        dem_raw = download_to(raw_dir/"dem"/"dem.tif", _dem)
    dem_raw = _resolve_raster_file(dem_raw, ("dem.tif","DEM.tif","dem.vrt"))


    # Diagnostics captions
    awc_sig = sniff_signature(str(awc_raw)) if awc_raw.exists() and awc_raw.is_file() else ("DIR" if awc_raw.exists() and awc_raw.is_dir() else "MISSING")
    dem_sig = sniff_signature(str(dem_raw)) if dem_raw.exists() and dem_raw.is_file() else ("DIR" if dem_raw.exists() and dem_raw.is_dir() else "MISSING")
    st.caption(f"AWC path: {awc_raw} • size={(awc_raw.stat().st_size if awc_raw.exists() and awc_raw.is_file() else 'NA')} • sig={awc_sig}")
    if awc_sig == "DIR":
        st.caption(f"AWC directory contents: { _dir_hint(awc_raw.parent if awc_raw.suffix else awc_raw) }")
    st.caption(f"DEM path: {dem_raw} • size={(dem_raw.stat().st_size if dem_raw.exists() and dem_raw.is_file() else 'NA')} • sig={dem_sig}")
    if dem_sig == "DIR":
        st.caption(f"DEM directory contents: { _dir_hint(dem_raw.parent if dem_raw.suffix else dem_raw) }")

    # Validate files
    ok_awc, msg_awc = (False, "Path is a directory; please point to a .tif/.vrt") if awc_raw.is_dir() else validate_geotiff(str(awc_raw))
    ok_dem, msg_dem = (False, "Path is a directory; please point to a .tif/.vrt") if dem_raw.is_dir() else validate_geotiff(str(dem_raw))
    if not ok_awc:
        st.error(f"AWC invalid: {msg_awc}. Ensure the link points to the actual .tif and is shared 'Anyone with link'.")
        st.stop()
    if not ok_dem:
        st.error(f"DEM invalid: {msg_dem}.")
        st.stop()

    S["awc_path"], S["dem_path"]=str(awc_raw), str(dem_raw)
    S["awc_loaded"], S["dem_loaded"] = True, True

    # Overlap diagnostics
    with rio.open(str(dem_raw)) as ds:
        dem_bounds = ds.bounds; dem_crs = ds.crs
    aoi = gpd.read_file(aoi_path)
    aoi_r = aoi.to_crs(dem_crs) if str(aoi.crs) != str(dem_crs) else aoi
    overlap = box(*dem_bounds).intersects(aoi_r.geometry.unary_union)
    st.caption(f"DEM bounds: {tuple(round(v,3) for v in dem_bounds)} | AOI overlaps DEM: {overlap}")

    # Clip & compute
    awc_da_clip = clip_to_aoi(awc_raw, aoi_path)

    try:
        dem_da_clip = clip_to_aoi(dem_raw, aoi_path)
        slope_da    = slope_from_dem(dem_raw)
        slope_da    = reproject_match(slope_da, dem_da_clip, resampling=Resampling.bilinear)
        dem_ok = True
    except Exception as e:
        st.warning(f"DEM could not be clipped to AOI ({e}). Falling back to neutral slope (0.5).")
        slope_da = xr.full_like(awc_da_clip, 0.5)
        dem_ok = False
        S["dem_loaded"] = False

    w_awc=float(st.session_state["awc_weight"]); w_slp=float(st.session_state["slope_weight"]); w_lc=float(st.session_state.get("lc_weight", cfg["weights"]["water"]["landcover"]))
    score = infiltration_score(awc_da_clip, slope_da, None, w_awc, w_slp, w_lc)

    # Ensure georeferencing for save
    def _ensure_like(ref: xr.DataArray, out, name="InfiltrationScore"):
        if isinstance(out, xr.DataArray) and tuple(out.dims)==tuple(ref.dims) and all(d in out.coords for d in ref.dims):
            da = out.copy(); da.name=name
        else:
            da = xr.DataArray(np.asarray(out), coords={ref.dims[0]: ref.coords[ref.dims[0]], ref.dims[1]: ref.coords[ref.dims[1]]}, dims=ref.dims, attrs=ref.attrs, name=name)
        if getattr(da.rio, "crs", None) is None and getattr(ref.rio, "crs", None) is not None:
            da = da.rio.write_crs(ref.rio.crs)
        try: _ = da.rio.transform(recalc=True)
        except Exception: da = da.rio.write_transform(ref.rio.transform(recalc=True))
        return da

    score = _ensure_like(awc_da_clip, score, name="InfiltrationScore")
    score = score.where(np.isfinite(score), np.nan)
    score.rio.write_nodata(np.nan, inplace=True)

    out_dir = Path(cfg["paths"]["processed"]) / "hydrology"; out_dir.mkdir(parents=True, exist_ok=True)
    out_tif = out_dir / "InfiltrationScore.tif"
    score.rio.to_raster(out_tif, compress="LZW")
    S["infiltration_path"]=str(out_tif)

    if dem_ok: st.success("✅ Analysis completed. Open the **Mapa** page to view the Infiltration layer.")
    else:      st.success("✅ Analysis completed (with neutral slope fallback). Open the **Mapa** page to view the Infiltration layer.")
