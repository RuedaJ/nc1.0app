"""
Map page: shows base map, AOI overlay, and tiled infiltration raster.
Includes a click-inspector that samples multiple rasters.
"""

from __future__ import annotations
from pathlib import Path
import streamlit as st
import folium
from streamlit_folium import st_folium

from app.tiles import raster_tile_layer
from etl.io import read_vector
from etl import hydrology as hydro


def repo_root_from_file(__file__: str) -> Path:
    return Path(__file__).resolve().parents[2]


_REPO_ROOT = repo_root_from_file(__file__)
DATA_DIR = _REPO_ROOT / "data"
PROC_DIR = DATA_DIR / "processed"
HYDRO_DIR = PROC_DIR / "hydrology"
AOI_DIR = PROC_DIR / "aoi"

st.set_page_config(page_title="Mapa", layout="wide")
st.title("🗺️ Mapa")


# ---------------------------
# State & inputs
# ---------------------------

state = st.session_state
lat = state.get("lat", 40.4168)
lon = state.get("lon", -3.7038)

paths = state.get("paths", {
    "awc": str(HYDRO_DIR / "AWC_clip.tif"),
    "dem": str(HYDRO_DIR / "DEM_clip.tif"),
    "slope": str(HYDRO_DIR / "Slope_deg.tif"),
    "infiltration": str(HYDRO_DIR / "InfiltrationScore.tif"),
})
# Persist back to session
state["paths"] = paths


# ---------------------------
# Build map
# ---------------------------

m = folium.Map(location=[lat, lon], zoom_start=12, tiles=None)
folium.TileLayer("CartoDB positron", name="Base").add_to(m)

# Tiled infiltration raster (if present)
infil_tif = Path(paths.get("infiltration", ""))
if infil_tif.exists():
    try:
        tile, _client = raster_tile_layer(infil_tif)
        tile.add_to(m)
    except Exception as e:
        st.warning(f"No se pudo usar tiles locales: {e}")

# AOI overlay (if present). Expect a GeoPackage: data/processed/aoi/aoi.gpkg
aoi_gpkg = AOI_DIR / "aoi.gpkg"
if aoi_gpkg.exists():
    try:
        gdf = read_vector(aoi_gpkg)
        if getattr(gdf, "crs", None) and gdf.crs.to_epsg() != 4326:
            gdf = gdf.to_crs(4326)
        folium.GeoJson(gdf.__geo_interface__, name="AOI").add_to(m)
    except Exception as e:
        st.warning(f"No se pudo mostrar AOI: {e}")

folium.LayerControl().add_to(m)

# Render map and capture interactions
st.write("Haz clic en el mapa para inspeccionar valores.")
event = st_folium(m, width=None, height=700, returned_objects=["last_clicked"])

# Inspector sampling
if event and isinstance(event, dict) and event.get("last_clicked"):
    ll = event["last_clicked"]
    click_lon, click_lat = float(ll["lng"]), float(ll["lat"])
    layers = {
        "AWC": paths.get("awc", ""),
        "Slope": paths.get("slope", ""),
        "Infiltration": paths.get("infiltration", ""),
    }
    vals = hydro.sample_rasters_at(click_lon, click_lat, {k.lower(): v for k, v in layers.items()})
    c1, c2, c3 = st.columns(3)
    c1.metric("AWC", f"{vals.get('awc', float('nan')):.3f}")
    c2.metric("Slope (°)", f"{vals.get('slope', float('nan')):.2f}")
    c3.metric("Infiltration", f"{vals.get('infiltration', float('nan')):.3f}")
