import streamlit as st

DEFAULTS = {
    "coords": None,
    "radius_km": 5,
    "aoi_geojson": None,
    "awc_link": None,
    "dem_link": None,
    "lc_link": None,
    "awc_path": None,
    "dem_path": None,
    "lc_path": None,
    "infiltration_path": None,
    "awc_loaded": False,
    "dem_loaded": False,
    "api_ok": None,
    "errors": [],
}

def get_state():
    for k, v in DEFAULTS.items():
        st.session_state.setdefault(k, v)
    return st.session_state

def reset():
    keys = list(DEFAULTS.keys())
    for k in keys:
        if k in st.session_state:
            st.session_state.pop(k)
