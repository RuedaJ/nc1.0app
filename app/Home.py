import sys, os
# Ensure repo root is on sys.path to import 'etl' package
_REPO_ROOT = os.path.dirname(os.path.dirname(__file__))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import streamlit as st
import yaml
from pathlib import Path
from state import get_state, reset
from etl.validators import valid_lat, valid_lon, looks_like_drive_url
from etl.geometry import aoi_from_latlon, save_geojson
from etl.drive import local_target_for

st.set_page_config(page_title="sustai-geo-app", layout="wide")

cfg = yaml.safe_load(Path("configs/default.yaml").read_text())
S = get_state()

st.sidebar.markdown("### 🌍 sustai-geo-app")
lang = st.sidebar.selectbox("Language/Idioma", ["EN","ES"], index=0 if cfg["app"]["language_default"]=="EN" else 1, key="lang_selector")

st.sidebar.markdown("---")
st.sidebar.markdown("**Site Location**")
lat = st.sidebar.text_input("Latitude", key="lat_input", placeholder="e.g., 41.6488")
lon = st.sidebar.text_input("Longitude", key="lon_input", placeholder="e.g., -0.8891")
use_map_click = st.sidebar.button("📍 Use map click", key="map_click_btn")

st.sidebar.markdown("**Analysis Area**")
radius_km = st.sidebar.slider("Analysis Radius (km)", min_value=1, max_value=cfg["app"]["max_radius_km"], value=S.get("radius_km", cfg["app"]["default_radius_km"]), key="radius_slider")

st.sidebar.markdown("**Data Sources (Google Drive)**")
awc_link = st.sidebar.text_input("AWC Raster Link *", value=S.get("awc_link") or "", key="awc_link_input", help="Google Drive shareable link")
dem_link = st.sidebar.text_input("DEM (PNOA) Link *", value=S.get("dem_link") or "", key="dem_link_input", help="Google Drive shareable link")
lc_link  = st.sidebar.text_input("Land Cover (optional)", value=S.get("lc_link") or "", key="lc_link_input")

with st.sidebar.expander("Advanced Settings"):
    st.slider("AWC Weight", 0.0, 1.0, cfg["weights"]["water"]["awc"], 0.05, key="awc_weight")
    st.slider("Slope Weight", 0.0, 1.0, cfg["weights"]["water"]["slope"], 0.05, key="slope_weight")
    st.slider("Land Cover Weight", 0.0, 1.0, cfg["weights"]["water"]["landcover"], 0.05, key="lc_weight")

run = st.sidebar.button("🚀 Run Analysis", use_container_width=True, key="run_analysis")
clear_cache = st.sidebar.button("🗑️ Clear Cache", key="clear_cache")

if clear_cache:
    reset()
    st.sidebar.success("Cache cleared and session reset.")

status = []
status.append("AWC ✓" if S.get("awc_loaded") else "AWC ▫")
status.append("DEM ✓" if S.get("dem_loaded") else "DEM ▫")
status.append("APIs ✓" if S.get("api_ok") else "APIs ▫")
st.sidebar.markdown("---")
st.sidebar.caption("Status: " + "  ".join(status))

st.title("Welcome 👋")
st.write("Use the tabs to navigate: **Mapa**, **Zonificación**, **Dashboard**, **Informe**.")
st.info("Enter coordinates and Google Drive links in the sidebar, then click **Run Analysis**. Data will be cached locally.")

if run:
    errors = []
    if not valid_lat(lat) or not valid_lon(lon):
        errors.append("Invalid coordinates. Please enter valid latitude and longitude.")
    if cfg["inputs"]["require_drive_links"]:
        _awc = st.session_state.get("awc_link_input") or ""
        _dem = st.session_state.get("dem_link_input") or ""
        if not looks_like_drive_url(_awc):
            errors.append("AWC: Invalid Google Drive link.")
        if not looks_like_drive_url(_dem):
            errors.append("DEM: Invalid Google Drive link.")
    if errors:
        for e in errors:
            st.error(e)
        S["errors"] = errors
    else:
        S["coords"] = (float(lat), float(lon))
        S["radius_km"] = int(radius_km)
        S["awc_link"] = st.session_state.get("awc_link_input") or ""
        S["dem_link"] = st.session_state.get("dem_link_input") or ""
        S["lc_link"] = (st.session_state.get("lc_link_input") or None)

        aoi_gdf = aoi_from_latlon(float(lat), float(lon), float(radius_km), cfg["app"]["crs"])
        aoi_path = Path(cfg["paths"]["processed"]) / "site_aoi" / f"site_buffer_{radius_km}km.geojson"
        save_geojson(aoi_gdf, aoi_path)
        S["aoi_geojson"] = str(aoi_path)

        raw_dir = Path(cfg["paths"]["raw"])
        S["awc_path"] = str(local_target_for(awc_link, raw_dir / "awc", "awc.tif"))
        S["dem_path"] = str(local_target_for(dem_link, raw_dir / "dem", "dem.tif"))
        if lc_link:
            S["lc_path"] = str(local_target_for(lc_link, raw_dir / "lc", "lc.gpkg"))

        st.success("Inputs validated and AOI created. Proceed to Mapa to visualize the AOI.")

st.markdown("**Design Spec Snapshot**")
st.write("- In MVP we focus on the **Water** module using AWC + DEM + (optional) land cover.")
st.write("- Section 0 groundwork is active: shared state, AOI builder, input validation, and caching policy hooks.")
