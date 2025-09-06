import sys, os
_APP_DIR = os.path.dirname(os.path.dirname(__file__))
if _APP_DIR not in sys.path:
    sys.path.insert(0, _APP_DIR)
from state import get_state

import streamlit as st
import yaml
from pathlib import Path
import leafmap.foliumap as leafmap


st.set_page_config(page_title="Mapa", layout="wide")

cfg = yaml.safe_load(Path("configs/default.yaml").read_text())
S = get_state()

st.title("Mapa")

col_map, col_layers = st.columns([4,1])

with col_map:
    m = leafmap.Map(center=cfg["app"]["map_center"], zoom=6, draw_control=False, locate_control=True)

    aoi_path = S.get("aoi_geojson")
    if aoi_path and Path(aoi_path).exists():
        m.add_geojson(aoi_path, layer_name="AOI (buffer)")
        try:
            m.zoom_to_bounds(aoi_path)
        except Exception:
            pass

    if S.get("coords"):
        lat, lon = S["coords"]
        m.add_marker(location=[lat, lon], popup="Site", draggable=False)

    infil_path = S.get("infiltration_path")
    if infil_path and Path(infil_path).exists():
        m.add_raster(str(infil_path), palette="Blues", layer_name="Infiltration")

    st.write(m.to_streamlit(height=600))

with col_layers:
    st.subheader("Layers")
    st.checkbox("Infiltration Score", value=bool(S.get("infiltration_path")), key="infiltration_toggle")
    st.slider("Opacity", 0.0, 1.0, 0.7, 0.05, key="infiltration_opacity")
    st.checkbox("Runoff Mask", value=False, key="runoff_toggle")
    st.checkbox("Flood T=500", value=False, key="flood_toggle")
    st.checkbox("Subcatchment", value=True, key="subcatchment_toggle")
    st.markdown("---")
    st.checkbox("Natura 2000", value=False, key="natura_toggle")
    st.checkbox("Slope Classes", value=False, key="slope_class_toggle")

st.markdown("#### Inspector")
st.caption("Click the map (in the final build) to see values: Infiltration, AWC, slope, land cover.")
