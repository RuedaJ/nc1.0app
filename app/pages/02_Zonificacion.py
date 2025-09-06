import sys, os
_APP_DIR = os.path.dirname(os.path.dirname(__file__))
if _APP_DIR not in sys.path:
    sys.path.insert(0, _APP_DIR)
from state import get_state


import streamlit as st

st.set_page_config(page_title="Zonificación", layout="wide")
S = get_state()

st.title("Zonificación")

c1, c2, c3 = st.columns(3)
with c1:
    method = st.radio("Classification Method", ["Natural Breaks", "Quantiles"], index=0, key="class_method")
with c2:
    excl_flood = st.checkbox("Exclude Floodplains", value=True, key="exclude_flood")
with c3:
    excl_prot  = st.checkbox("Exclude Protected Areas", value=False, key="exclude_protected")

if not S.get("infiltration_path"):
    st.info("Run analysis to compute the Infiltration Score, then zoning will be available here.")

c4, c5 = st.columns(2)
with c4:
    st.download_button("📥 GeoJSON", b"", file_name="priority_zones.geojson", disabled=True, key="export_geojson")
with c5:
    st.download_button("📥 GeoPackage", b"", file_name="priority_zones.gpkg", disabled=True, key="export_gpkg")

m1, m2, m3, m4 = st.columns(4)
m1.metric("Priority A", "0 ha")
m2.metric("Priority B", "0 ha")
m3.metric("Priority C", "0 ha")
m4.write("**Recommended Actions**\n\n• Green infrastructure\n\n• Rain gardens\n\n• Permeable surfaces")
