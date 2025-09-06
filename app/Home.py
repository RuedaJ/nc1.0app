import sys, os
# Ensure repo root is on sys.path to import local 'etl' package and 'state'
_REPO_ROOT = os.path.dirname(os.path.dirname(__file__))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import streamlit as st
import yaml
from pathlib import Path

from state import get_state, reset

# Validators & Drive helpers
from etl.validators import (
    valid_lat, valid_lon, looks_like_drive_url,
    validate_geotiff, sniff_signature
)
from etl.drive import extract_id, download_drive_file, download_to
from etl.geometry import aoi_from_latlon, save_geojson

# Hydrology ops
from etl.hydrology import (
    clip_to_aoi, slope_from_dem, reproject_match, infiltration_score
)
from rasterio.enums import Resampling

st.set_page_config(page_title="sustai-geo-app", layout="wide")

# ───────────────────── Config & State ─────────────────────
cfg = yaml.safe_load(Path("configs/default.yaml").read_text())
S = get_state()

# ───────────────────── Sidebar ─────────────────────
st.sidebar.markdown("### 🌍 sustai-geo-app")
lang = st.sidebar.selectbox(
    "Language/Idioma", ["EN", "ES"],
    index=0 if cfg["app"]["language_default"] == "EN" else 1,
    key="lang_selector"
)

st.sidebar.markdown("---")
st.sidebar.markdown("**Site Location**")
lat = st.sidebar.text_input("Latitude", key="lat_input", placeholder="e.g., 41.6488")
lon = st.sidebar.text_input("Longitude", key="lon_input", placeholder="e.g., -0.8891")
st.sidebar.button("📍 Use map click", key="map_click_btn")

st.sidebar.markdown("**Analysis Area**")
radius_km = st.sidebar.slider(
    "Analysis Radius (km)",
    min_value=1,
    max_value=cfg["app"]["max_radius_km"],
    value=S.get("radius_km", cfg["app"]["default_radius_km"]),
    key="radius_slider"
)

st.sidebar.markdown("**Data Sources (Google Drive or local path)**")
st.sidebar.text_input("AWC Raster Link *", value=S.get("awc_link") or "",
                      key="awc_link_input", help="Google Drive shareable link OR a local path to a .tif")
st.sidebar.text_input("DEM (PNOA) Link *", value=S.get("dem_link") or "",
                      key="dem_link_input", help="Google Drive shareable link OR a local path to a .tif")
st.sidebar.text_input("Land Cover (optional)", value=S.get("lc_link") or "",
                      key="lc_link_input")

with st.sidebar.expander("Advanced Settings"):
    st.slider("AWC Weight", 0.0, 1.0, cfg["weights"]["water"]["awc"], 0.05, key="awc_weight")
    st.slider("Slope Weight", 0.0, 1.0, cfg["weights"]["water"]["slope"], 0.05, key="slope_weight")
    # IMPORTANT: we consistently use this key later
    st.slider("Land Cover Weight", 0.0, 1.0, cfg["weights"]["water"]["landcover"], 0.05, key="lc_weight")

run = st.sidebar.button("🚀 Run Analysis", use_container_width=True, key="run_analysis")
if st.sidebar.button("🗑️ Clear Cache", key="clear_cache"):
    reset()
    st.sidebar.success("Cache cleared and session reset.")

status = []
status.append("AWC ✓" if S.get("awc_loaded") else "AWC ▫")
status.append("DEM ✓" if S.get("dem_loaded") else "DEM ▫")
status.append("APIs ✓" if S.get("api_ok") else "APIs ▫")
st.sidebar.markdown("---")
st.sidebar.caption("Status: " + "  ".join(status))

# ───────────────────── Main ─────────────────────
st.title("Welcome 👋")
st.write("Use the tabs to navigate: **Mapa**, **Zonificación**, **Dashboard**, **Informe**.")
st.info("Enter coordinates and Google Drive links (or local paths) in the sidebar, then click **Run Analysis**. Data will be cached locally.")

# ───────────────────── Run Analysis ─────────────────────
if run:
    errors = []

    # 1) Validate inputs
    if not valid_lat(lat) or not valid_lon(lon):
        errors.append("Invalid coordinates. Please enter valid latitude and longitude.")
    _awc = st.session_state.get("awc_link_input") or ""
    _dem = st.session_state.get("dem_link_input") or ""
    _lc  = st.session_state.get("lc_link_input") or ""

    if cfg["inputs"]["require_drive_links"]:
        if not (Path(_awc).exists() or looks_like_drive_url(_awc)):
            errors.append("AWC: Invalid Google Drive link.")
        if not (Path(_dem).exists() or looks_like_drive_url(_dem)):
            errors.append("DEM: Invalid Google Drive link.")

    if errors:
        for e in errors:
            st.error(e)
        S["errors"] = errors
        st.stop()

    # 2) Save inputs to state
    S["coords"] = (float(lat), float(lon))
    S["radius_km"] = int(radius_km)
    S["awc_link"] = _awc
    S["dem_link"] = _dem
    S["lc_link"]  = (_lc or None)

    # 3) Build AOI
    aoi_gdf = aoi_from_latlon(float(lat), float(lon), float(radius_km), cfg["app"]["crs"])
    aoi_path = Path(cfg["paths"]["processed"]) / "site_aoi" / f"site_buffer_{radius_km}km.geojson"
    save_geojson(aoi_gdf, aoi_path)
    S["aoi_geojson"] = str(aoi_path)

    # 4) Download/resolve AWC & DEM (robust)
    raw_dir = Path(cfg["paths"]["raw"])

    # AWC
    if Path(_awc).exists():
        awc_raw = Path(_awc)  # local path fallback
    elif extract_id(_awc):  # Google Drive link
        awc_raw = download_drive_file(_awc, raw_dir / "awc")  # preserves Drive filename
    else:
        awc_raw = download_to(raw_dir / "awc" / "awc.tif", _awc)

    # DEM
    if Path(_dem).exists():
        dem_raw = Path(_dem)
    elif extract_id(_dem):
        dem_raw = download_drive_file(_dem, raw_dir / "dem")
    else:
        dem_raw = download_to(raw_dir / "dem" / "dem.tif", _dem)

    # 5) Validate rasters before any processing
    ok_awc, msg_awc = validate_geotiff(str(awc_raw))
    ok_dem, msg_dem = validate_geotiff(str(dem_raw))
    st.caption(f"AWC path: {awc_raw} • size={awc_raw.stat().st_size if Path(awc_raw).exists() else 'NA'} • sig={sniff_signature(str(awc_raw))}")
    st.caption(f"DEM path: {dem_raw} • size={dem_raw.stat().st_size if Path(dem_raw).exists() else 'NA'} • sig={sniff_signature(str(dem_raw))}")
    if not ok_awc:
        st.error(f"AWC invalid: {msg_awc}. Ensure the link points to the actual .tif and is shared 'Anyone with link'.")
        st.stop()
    if not ok_dem:
        st.error(f"DEM invalid: {msg_dem}.")
        st.stop()

    S["awc_path"], S["dem_path"] = str(awc_raw), str(dem_raw)
    S["awc_loaded"], S["dem_loaded"] = True, True

    # 6) Clip to AOI & compute slope
    awc_da_clip = clip_to_aoi(awc_raw, aoi_path)       # AWC in AOI
    dem_da_clip = clip_to_aoi(dem_raw, aoi_path)       # DEM in AOI
    slope_da    = slope_from_dem(dem_raw)              # slope on full DEM grid
    slope_da    = reproject_match(slope_da, dem_da_clip, resampling=Resampling.bilinear)

    # 7) Compute Infiltration Score (LC optional; fallback handled in function when None)
    lc_da = None
    w_awc = float(st.session_state["awc_weight"])
    w_slp = float(st.session_state["slope_weight"])
    w_lc  = float(st.session_state.get("lc_weight", cfg["weights"]["water"]["landcover"]))

    score = infiltration_score(awc_da_clip, slope_da, lc_da, w_awc, w_slp, w_lc)

    # 8) Save raster & set state
    out_dir = Path(cfg["paths"]["processed"]) / "hydrology"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tif = out_dir / "InfiltrationScore.tif"
    score.rio.to_raster(out_tif, compress="LZW")
    S["infiltration_path"] = str(out_tif)

    st.success("✅ Analysis completed. Open the **Mapa** page to view the Infiltration layer.")

# ───────────────────── Footer ─────────────────────
st.markdown("**Design Spec Snapshot**")
st.write("- In MVP we focus on the **Water** module using AWC + DEM + (optional) land cover.")
st.write("- Section 0 groundwork is active: shared state, AOI builder, input validation, and caching policy hooks.")
