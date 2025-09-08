"""
Dashboard page: quick KPIs and system health.
"""

from __future__ import annotations
from pathlib import Path
import json
import streamlit as st
import psutil


def repo_root_from_file(__file__: str) -> Path:
    return Path(__file__).resolve().parents[2]


_REPO_ROOT = repo_root_from_file(__file__)
DATA_DIR = _REPO_ROOT / "data"
PROC_DIR = DATA_DIR / "processed"
HYDRO_DIR = PROC_DIR / "hydrology"

st.set_page_config(page_title="Dashboard", layout="wide")
st.title("📊 Dashboard")


# ---------------------------
# System health
# ---------------------------

vm = psutil.virtual_memory()
c1, c2, c3 = st.columns(3)
with c1:
    st.metric("RAM usada", f"{(vm.total - vm.available)/1e9:.1f} / {vm.total/1e9:.1f} GB")
with c2:
    st.metric("CPU", f"{psutil.cpu_percent(interval=0.3):.0f}%")
with c3:
    st.metric("Procesos", f"{len(psutil.pids())}")


# ---------------------------
# KPIs from session or compute lightweight stats
# ---------------------------

state = st.session_state
paths = state.get("paths", {})
run_stats = state.get("run_stats", {})

st.subheader("KPIs (AOI)")
k_cols = st.columns(4)

def _calc_mean_from_tif(path: str | Path) -> float | None:
    try:
        import rasterio
        import numpy as np
        with rasterio.open(path) as src:
            arr = src.read(1, masked=True)
        if arr.size == 0:
            return None
        return float(np.nanmean(arr))
    except Exception:
        return None

awc_mean = _calc_mean_from_tif(paths.get("awc")) if paths.get("awc") else None
slope_mean = _calc_mean_from_tif(paths.get("slope")) if paths.get("slope") else None
infil_mean = _calc_mean_from_tif(paths.get("infiltration")) if paths.get("infiltration") else None

with k_cols[0]:
    st.metric("AWC (media)", f"{awc_mean:.3f}" if awc_mean is not None else "—")
with k_cols[1]:
    st.metric("Pendiente (° media)", f"{slope_mean:.2f}" if slope_mean is not None else "—")
with k_cols[2]:
    st.metric("Infiltración (media)", f"{infil_mean:.3f}" if infil_mean is not None else "—")
with k_cols[3]:
    if run_stats:
        st.metric("Duración (s)", run_stats.get("duration_s", "—"))
    else:
        st.metric("Duración (s)", "—")

st.divider()

st.subheader("Archivos de salida")
table_data = []
for label, key in [("AWC (clip)", "awc"), ("DEM (clip)", "dem"), ("Slope", "slope"), ("Infiltration", "infiltration")]:
    p = paths.get(key)
    table_data.append({
        "Capa": label,
        "Ruta": str(p) if p else "—",
        "Existe": Path(p).exists() if p else False,
        "Tamaño (MB)": round(Path(p).stat().st_size/1e6, 2) if p and Path(p).exists() else 0,
    })

st.dataframe(table_data, use_container_width=True)
