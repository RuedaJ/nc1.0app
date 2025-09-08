"""
Report page: compose a simple report, show branded map image, export downloads.
"""

from __future__ import annotations
# === PATH BOOTSTRAP ===
import sys
from pathlib import Path
_THIS = Path(__file__).resolve()
_REPO_ROOT = _THIS.parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))
# ======================

import io
import datetime as dt
import streamlit as st
from app.imaging import add_footer_badge, thumbnail

DATA_DIR = _REPO_ROOT / "data"
PROC_DIR = DATA_DIR / "processed"
HYDRO_DIR = PROC_DIR / "hydrology"
REPORT_DIR = DATA_DIR / "reports"
REPORT_DIR.mkdir(parents=True, exist_ok=True)

st.set_page_config(page_title="Informe", layout="wide")
st.title("🧾 Informe")

# ---------------------------
# Load state / inputs
# ---------------------------
state = st.session_state
paths = state.get("paths", {})  # keys: awc, dem, slope, infiltration
run_stats = state.get("run_stats", {})

st.subheader("Resumen")
cols = st.columns(2)
with cols[0]:
    st.write("- **Fecha:**", dt.datetime.now().strftime("%Y-%m-%d %H:%M"))
    st.write("- **Resultados:**", "InfiltrationScore.tif" if paths.get("infiltration") else "—")
with cols[1]:
    if run_stats:
        st.write(
            f"- **Duración:** {run_stats.get('duration_s', '—')} s  \n"
            f"- **RAM pico:** {run_stats.get('peak_gb', '—')} GB"
        )

st.divider()

st.subheader("Mapa (brandeado)")
# If you have a map snapshot saved elsewhere, use it; otherwise a placeholder is created.
placeholder_img = REPORT_DIR / "map_placeholder.png"
if not placeholder_img.exists():
    try:
        from PIL import Image, ImageDraw
        im = Image.new("RGB", (1200, 700), (245, 247, 250))
        d = ImageDraw.Draw(im)
        d.text((20, 20), "Map placeholder — add your snapshot export here", fill=(60, 60, 60))
        im.save(placeholder_img)
    except Exception:
        pass

base_img = placeholder_img if placeholder_img.exists() else None
if base_img:
    branded = add_footer_badge(base_img, REPORT_DIR / "map_branded.jpg", text="sustai-geo-app")
    thumb = thumbnail(branded, REPORT_DIR / "map_thumb.jpg", max_side=480)
    st.image(str(branded), caption="Mapa (con banda de marca)")
    with open(branded, "rb") as fh:
        st.download_button("Descargar mapa (JPG)", data=fh.read(), file_name=Path(branded).name)

st.divider()

st.subheader("Fuentes y atribuciones")
with st.expander("Ver detalles"):
    st.markdown(
        """
- **AWC**: Google Drive (enlace proporcionado por el usuario)
- **DEM**: Google Drive (enlace proporcionado por el usuario)
- **Procesos**: recorte AOI, reproyección, pendiente (DEM), puntuación de infiltración.
- **Attribution**: PNOA / IGN España (si aplicable), demás fuentes según origen de datos.
        """
    )

st.divider()

st.subheader("Descargas")
dl_cols = st.columns(4)
for i, (label, key) in enumerate([("AWC (clip)", "awc"), ("DEM (clip)", "dem"), ("Slope", "slope"), ("Infiltration", "infiltration")]):
    with dl_cols[i]:
        p = paths.get(key)
        if p and Path(p).exists():
            with open(p, "rb") as fh:
                st.download_button(f"Descargar {label}", data=fh.read(), file_name=Path(p).name)
        else:
            st.caption(f"{label}: —")

st.info("Para un informe completo, renderiza Markdown/PDF con WeasyPrint usando plantillas Jinja2 (próximo paso).")
