"""
Mapa: base + AOI + Infiltration raster with opacity and inspector.
Uses localtileserver when available; falls back to ImageOverlay otherwise.
"""

from __future__ import annotations

# === PATH BOOTSTRAP (ensures project root is importable) ===
import sys
from pathlib import Path
_THIS = Path(__file__).resolve()
_REPO_ROOT = _THIS.parents[2]  # project root (../../ from app/pages/)
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))
# ==========================================================

import io
import numpy as np
import streamlit as st
import folium
from streamlit_folium import st_folium

# Try to use tiled rasters if helper exists; otherwise, no crash.
try:
    from app.tiles import raster_tile_layer
except Exception:
    raster_tile_layer = None  # fallback later

from etl.io import read_vector
from etl import hydrology as hydro

# Optional bits for the static ImageOverlay fallback
try:
    import rasterio
    from rasterio.warp import transform_bounds
except Exception:  # keep page working even if rasterio import fails
    rasterio = None
    transform_bounds = None

try:
    import branca
    from branca.colormap import linear
except Exception:
    branca = None
    linear = None


# ---------------------------
# Paths
# ---------------------------
DATA_DIR = _REPO_ROOT / "data"
PROC_DIR = DATA_DIR / "processed"
HYDRO_DIR = PROC_DIR / "hydrology"
AOI_DIR = PROC_DIR / "aoi"
AOI_GPKG = AOI_DIR / "aoi.gpkg"

DEFAULT_PATHS = {
    "awc": str(HYDRO_DIR / "AWC_clip.tif"),
    "dem": str(HYDRO_DIR / "DEM_clip.tif"),
    "slope": str(HYDRO_DIR / "Slope_deg.tif"),
    "infiltration": str(HYDRO_DIR / "InfiltrationScore.tif"),
}

# ---------------------------
# Page setup
# ---------------------------
st.set_page_config(page_title="Mapa", layout="wide")
st.title("🗺️ Mapa")

state = st.session_state
paths = state.get("paths", DEFAULT_PATHS.copy())
state["paths"] = paths  # persist
lat = state.get("lat", 40.4168)
lon = state.get("lon", -3.7038)

# ---------------------------
# Sidebar controls
# ---------------------------
with st.sidebar:
    st.header("Capas")
    basemap = st.selectbox(
        "Base",
        ["CartoDB positron", "OpenStreetMap", "Stamen Terrain", "Stamen Toner"],
        index=0,
    )

    show_infil = st.checkbox("Infiltration", value=True)
    infil_opacity = st.slider("Opacidad Infiltration", min_value=0.0, max_value=1.0, value=0.85, step=0.05)
    show_aoi = st.checkbox("AOI", value=True)

    st.divider()
    zoom_to_aoi = st.button("Zoom a AOI")

    st.divider()
    st.caption("Consejo: si no hay tiles locales, se usará una superposición estática.")

# ---------------------------
# Map build
# ---------------------------
m = folium.Map(location=[lat, lon], zoom_start=12, tiles=None)
folium.TileLayer(basemap, name="Base").add_to(m)

# ---- AOI layer
def _add_aoi(map_obj: folium.Map):
    if AOI_GPKG.exists():
        try:
            gdf = read_vector(AOI_GPKG)
            if getattr(gdf, "crs", None) and gdf.crs.to_epsg() != 4326:
                gdf = gdf.to_crs(4326)
            gj = folium.GeoJson(gdf.__geo_interface__, name="AOI", style_function=lambda x: {
                "color": "#1f77b4",
                "weight": 2,
                "fillOpacity": 0.05,
            })
            gj.add_to(map_obj)
            return gdf
        except Exception as e:
            st.warning(f"No se pudo mostrar AOI: {e}")
    else:
        st.caption("AOI: (no encontrado)")
    return None

gdf_aoi = None
if show_aoi:
    gdf_aoi = _add_aoi(m)

# ---- Infiltration raster layer
def _add_infiltration_tiles(map_obj: folium.Map, tif_path: Path, opacity: float) -> bool:
    """
    Try localtileserver tiles first; if not available, use a rasterio-based ImageOverlay fallback.
    Returns True if something was added, else False.
    """
    # Preferred: localtileserver
    if raster_tile_layer is not None:
        try:
            tile, _client = raster_tile_layer(tif_path)
            # Update opacity on the tile layer if supported
            try:
                tile.options["opacity"] = opacity
            except Exception:
                pass
            tile.add_to(map_obj)
            return True
        except Exception as e:
            st.warning(f"No se pudo usar tiles locales (localtileserver): {e}. Probando superposición estática...")
    # Fallback: static ImageOverlay using raster bounds -> WGS84
    if rasterio is None or transform_bounds is None:
        st.warning("Rasterio no disponible para la superposición estática.")
        return False
    try:
        with rasterio.open(tif_path) as src:
            # bounds in source CRS → WGS84 for Leaflet
            left, bottom, right, top = transform_bounds(src.crs, "EPSG:4326", *src.bounds, densify_pts=21)
            # Render a small, downsampled PNG in memory for overlay
            # (avoid huge images; this is a fallback)
            data = src.read(1, masked=True)
            # Normalize to [0, 255] for preview
            arr = np.array(data, dtype="float64")
            vmin = np.nanpercentile(arr, 2)
            vmax = np.nanpercentile(arr, 98)
            rng = (vmax - vmin) if vmax > vmin else 1.0
            norm = np.clip((arr - vmin) / rng, 0, 1)
            rgba = np.dstack([norm*255, norm*255, norm*255, np.where(np.isnan(arr), 0, opacity*255)]).astype("uint8")

            # Encode to PNG
            from PIL import Image
            buf = io.BytesIO()
            Image.fromarray(rgba, mode="RGBA").save(buf, format="PNG")
            buf.seek(0)

            folium.raster_layers.ImageOverlay(
                image=buf,
                bounds=[[bottom, left], [top, right]],
                name="Infiltration (static)",
                opacity=opacity,
                interactive=False,
                cross_origin=False,
                zindex=1,
            ).add_to(map_obj)
            return True
    except Exception as e:
        st.warning(f"Fallo en superposición estática: {e}")
        return False

infil_tif = Path(paths.get("infiltration", ""))
if show_infil and infil_tif.exists():
    _ = _add_infiltration_tiles(m, infil_tif, infil_opacity)
else:
    if show_infil:
        st.info("Infiltration: archivo no encontrado.")

# ---- Legend for infiltration (simple)
def _add_infiltration_legend(map_obj: folium.Map):
    if branca is None or linear is None:
        return
    cmap = linear.Greys_09.scale(0, 1.0)
    cmap.caption = "Infiltration (0–1)"
    map_obj.add_child(cmap)

if show_infil:
    _add_infiltration_legend(m)

# ---- Layer control
folium.LayerControl().add_to(m)

# ---- Zoom to AOI if requested
if zoom_to_aoi and gdf_aoi is not None and len(gdf_aoi) > 0:
    try:
        bounds = gdf_aoi.total_bounds  # minx, miny, maxx, maxy in 4326
        # Store it in session and the map will center after render via fit_bounds
        state["fit_bounds"] = [[bounds[1], bounds[0]], [bounds[3], bounds[2]]]
    except Exception as e:
        st.warning(f"No se pudo centrar a AOI: {e}")

# Render map and capture interactions
ret = st_folium(
    m,
    width=None,
    height=700,
    returned_objects=["last_clicked"],
    center=state.get("fit_bounds") is None,
    feature_group_to_add="",
)

# If we queued a fit_bounds, instruct the frontend. (st_folium supports center param but not direct fit post-render;
# this is a best-effort: the user can click zoom to AOI and we set the center on next rerun.)
if "fit_bounds" in state:
    # Clear it on the next run; Streamlit reruns will refresh.
    del state["fit_bounds"]

# ---------------------------
# Inspector sampling
# ---------------------------
if ret and isinstance(ret, dict) and ret.get("last_clicked"):
    ll = ret["last_clicked"]
    click_lon, click_lat = float(ll["lng"]), float(ll["lat"])
    layers = {
        "AWC": paths.get("awc", ""),
        "Slope": paths.get("slope", ""),
        "Infiltration": paths.get("infiltration", ""),
    }
    vals = hydro.sample_rasters_at(click_lon, click_lat, {k.lower(): v for k, v in layers.items()})

    c1, c2, c3 = st.columns(3)
    def _fmt(x, n=3):
        try:
            if np.isnan(x):
                return "—"
            return f"{x:.{n}f}"
        except Exception:
            return "—"

    c1.metric("AWC", _fmt(vals.get("awc", float("nan")), 3))
    c2.metric("Pendiente (°)", _fmt(vals.get("slope", float("nan")), 2))
    c3.metric("Infiltración", _fmt(vals.get("infiltration", float("nan")), 3))
