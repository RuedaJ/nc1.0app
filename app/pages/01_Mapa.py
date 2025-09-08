# app/pages/01_Mapa.py
import sys, os
from pathlib import Path

# Ensure repo root on path so we can import state & etl modules
_REPO_ROOT = os.path.dirname(os.path.dirname(__file__))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import streamlit as st
import geopandas as gpd
import numpy as np
import rasterio as rio
from rasterio.enums import Resampling
from rasterio.windows import Window
import rioxarray as rxr
from shapely.geometry import box
from pyproj import Transformer

import folium
from streamlit_folium import st_folium

# Color maps (avoid old deprecations)
try:
    from matplotlib import colormaps as _cmaps
    def _get_cmap(name: str):
        return _cmaps.get(name)
except Exception:
    import matplotlib.cm as _cmaps  # fallback
    def _get_cmap(name: str):
        return _cmaps.get_cmap(name)

from branca.colormap import LinearColormap

from state import get_state

st.set_page_config(page_title="Mapa — sustai-geo-app", layout="wide")
S = get_state()

st.title("Mapa")
st.caption("AOI, Infiltration layer, and click-to-inspect values.")

# ---- Guards ----
infil_path = S.get("infiltration_path")
aoi_geojson = S.get("aoi_geojson")
coords = S.get("coords")

if not infil_path or not Path(infil_path).exists():
    st.warning("No Infiltration layer found. Run the analysis on the Home page first.")
    st.stop()

# ---- Sidebar controls ----
with st.sidebar:
    st.markdown("### Layers")
    opacity = st.slider("Infiltration opacity", 0.0, 1.0, 0.75, 0.05)
    palette = st.selectbox("Palette", ["Blues", "viridis", "turbo", "magma"], index=0)
    stretch = st.selectbox("Stretch", ["min–max", "p2–p98"], index=1, help="How to normalize values for display.")
    show_aoi = st.checkbox("Show AOI", value=True)
    show_site = st.checkbox("Show site marker", value=True)
    show_overlay_outline = st.checkbox("Debug: show overlay outline", value=False)

# ---------- Helpers ----------
def add_aoi_to_folium_map(fm: folium.Map, aoi_path: str):
    try:
        gdf = gpd.read_file(aoi_path)
        gj = folium.GeoJson(
            gdf.to_crs(4326).__geo_interface__,
            name="AOI",
            style_function=lambda _:
                {"color": "#0066CC", "weight": 2, "fill": False},
        )
        gj.add_to(fm)
        # Keep last fit to the raster; so AOI fit happens only if raster is missing
        if "Infiltration" not in [x.layer_name if hasattr(x, "layer_name") else getattr(x, "overlay", None) for x in fm._children.values()]:
            minx, miny, maxx, maxy = gdf.to_crs(4326).total_bounds
            fm.fit_bounds([[miny, minx], [maxy, maxx]])
    except Exception as e:
        st.warning(f"Could not render AOI: {e}")

def add_site_marker_to_folium_map(fm: folium.Map, lat: float, lon: float):
    folium.CircleMarker(
        location=[lat, lon],
        radius=5,
        color="#DC3545",
        fill=True,
        fill_opacity=1,
        popup="Site",
    ).add_to(fm)

def _normalize_for_display(arr: np.ndarray, kind: str):
    finite = np.isfinite(arr)
    if not finite.any():
        st.caption("Infiltration display: no finite values; rendering blank.")
        return np.zeros_like(arr), 0.0, 0.0, 0.0

    if kind == "p2–p98":
        p2, p98 = np.nanpercentile(arr[finite], [2, 98])
        norm = np.clip((arr - p2) / max(p98 - p2, 1e-6), 0, 1)
        st.caption(f"Infiltration display: valid%={100*np.mean(finite):.1f} (p2={p2:.3f}, p98={p98:.3f})")
        return norm, float(p2), float(p98), 100*np.mean(finite)
    else:
        vmin, vmax = float(np.nanmin(arr[finite])), float(np.nanmax(arr[finite]))
        norm = (arr - vmin) / max(vmax - vmin, 1e-6)
        st.caption(f"Infiltration display: valid%={100*np.mean(finite):.1f}, min={vmin:.4f}, max={vmax:.4f}")
        return norm, vmin, vmax, 100*np.mean(finite)

def add_raster_to_folium(fm: folium.Map, tif_path: str, name: str, opacity: float, palette_name: str, stretch_kind: str, show_outline: bool):
    """
    Preferred: localtileserver (if installed) via leafmap.common.get_local_tile_layer
    Fallback: reproject -> (optional downsample) -> normalize -> colorize -> ImageOverlay with live opacity
    """
    # Try local tiles first (if available)
    try:
        from leafmap.common import get_local_tile_layer
        # Do NOT pass name= to avoid "multiple values" error; set name later if possible.
        tile = get_local_tile_layer(
            tif_path,
            palette=palette_name,
            vmin=0, vmax=1,
            opacity=float(opacity),
        )
        try:
            # Some objects are directly addable
            tile.add_to(fm)
        except Exception:
            # Others require add_child
            fm.add_child(tile)
        # Try to assign a friendly label for LayerControl
        try:
            setattr(tile, "layer_name", name)
        except Exception:
            pass
        folium.LayerControl(collapsed=False).add_to(fm)
        return
    except Exception as e:
        st.info(f"localtileserver not available ({e}). Using ImageOverlay fallback.")

    # Fallback pipeline
    da = rxr.open_rasterio(tif_path, masked=True).squeeze()
    da_4326 = da.rio.reproject("EPSG:4326", resampling=Resampling.bilinear)

    # Downsample for snappy UI
    max_px = 1024
    h = int(da_4326.sizes[da_4326.dims[0]])
    w = int(da_4326.sizes[da_4326.dims[1]])
    scale = min(1.0, max_px / max(h, w))
    if scale < 1.0:
        resx, resy = da_4326.rio.resolution()
        da_4326 = da_4326.rio.reproject(
            "EPSG:4326",
            resolution=(abs(resx) / scale, abs(resy) / scale),
            resampling=Resampling.bilinear,
        )

    arr = da_4326.values.astype("float32")
    norm, a, b, valid_pc = _normalize_for_display(arr, stretch_kind)

    cmap = _get_cmap(palette_name)
    rgba = (cmap(np.clip(norm, 0, 1)) * 255).astype("uint8")

    # Bounds come as (minx, miny, maxx, maxy) = (west, south, east, north)
    west, south, east, north = da_4326.rio.bounds()

    # Set origin to match array orientation: if first row represents north -> origin='upper'
    # In geographic coords increasing y is north; DataArray y usually decreases from top to bottom => origin='upper'
    origin = "upper"

    folium.raster_layers.ImageOverlay(
        image=rgba,                              # pass RGBA ndarray
        bounds=[[south, west], [north, east]],
        name=name,
        opacity=float(opacity),                  # bind slider value
        interactive=False,
        cross_origin=False,
        zindex=500,
        origin=origin,
    ).add_to(fm)

    # Optional: draw outline box for debug
    if show_outline:
        folium.Rectangle(
            bounds=[[south, west], [north, east]],
            color="#7B68EE",
            weight=1,
            fill=False,
            dash_array="4,4",
            tooltip="Infiltration extent",
        ).add_to(fm)

    # Fit to overlay bounds to ensure it's visible
    fm.fit_bounds([[south, west], [north, east]])
    folium.LayerControl(collapsed=False).add_to(fm)

    # Optional legend
    try:
        legend = LinearColormap(
            colors=[_get_cmap(palette_name)(x) for x in np.linspace(0, 1, 6)],
            vmin=0, vmax=1
        ).to_step(index=[0, 0.2, 0.4, 0.6, 0.8, 1])
        legend.caption = "Infiltration (0 → 1)"
        legend.add_to(fm)
    except Exception:
        pass

def sample_tiff_value(tif_path: str, lon: float, lat: float):
    try:
        with rio.open(tif_path) as ds:
            if ds.crs is None:
                return None
            transformer = Transformer.from_crs("EPSG:4326", ds.crs, always_xy=True)
            x, y = transformer.transform(lon, lat)
            val = list(ds.sample([(x, y)]))[0][0]
            if isinstance(val, np.ndarray):
                val = float(val[0])
            return None if (val is None or np.isnan(val)) else float(val)
    except Exception:
        return None

def sample_slope_from_dem(dem_path: str, lon: float, lat: float):
    """Approximate normalized slope [0,1] from a 3x3 window at click location."""
    try:
        if not dem_path or not Path(dem_path).exists():
            return None
        with rio.open(dem_path) as ds:
            if ds.crs is None:
                return None
            transformer = Transformer.from_crs("EPSG:4326", ds.crs, always_xy=True)
            x, y = transformer.transform(lon, lat)
            r, c = ds.index(x, y)
            win = Window(
                col_off=max(c - 1, 0),
                row_off=max(r - 1, 0),
                width=3 if c + 2 < ds.width else min(3, ds.width - max(c - 1, 0)),
                height=3 if r + 2 < ds.height else min(3, ds.height - max(r - 1, 0)),
            )
            arr = ds.read(1, window=win, boundless=True, fill_value=np.nan).astype("float32")
            if np.isnan(arr).all():
                return None
            a, _, _, _, e, _, _, _, _ = ds.transform
            gy, gx = np.gradient(arr)
            with np.errstate(invalid="ignore", divide="ignore"):
                slope_rad = np.arctan(np.sqrt((gx * (1.0 / abs(a))) ** 2 + (gy * (1.0 / abs(e))) ** 2))
            slope_deg = np.degrees(slope_rad)
            if np.isnan(slope_deg).all():
                return None
            smin, smax = np.nanmin(slope_deg), np.nanmax(slope_deg)
            if smax > smin:
                s_norm = (smax - slope_deg) / (smax - smin + 1e-6)
            else:
                s_norm = np.ones_like(slope_deg)
            return float(np.nanmean(s_norm))
    except Exception:
        return None

# ---------- Render ----------
base_loc = (coords[0], coords[1]) if coords else (41.65, -0.89)
fm = folium.Map(location=base_loc, zoom_start=12, control_scale=True, tiles="CartoDB positron")

# Raster first (so its fit_bounds wins)
add_raster_to_folium(fm, str(infil_path), "Infiltration", opacity, palette, stretch, show_overlay_outline)

# AOI & site
if aoi_geojson and Path(aoi_geojson).exists() and show_aoi:
    add_aoi_to_folium_map(fm, aoi_geojson)
if show_site and coords:
    add_site_marker_to_folium_map(fm, coords[0], coords[1])

# Render & click inspector with streamlit-folium
st_data = st_folium(fm, height=640, use_container_width=True, returned_objects=["last_clicked"])
last = st_data.get("last_clicked")
with st.expander("Inspector", expanded=bool(last)):
    if last:
        lat_click = last["lat"]; lon_click = last["lng"]
        infil_val = sample_tiff_value(infil_path, lon_click, lat_click)
        awc_path = S.get("awc_path")
        awc_val = sample_tiff_value(awc_path, lon_click, lat_click) if awc_path and Path(awc_path).exists() else None
        dem_path = S.get("dem_path")
        slope_val = sample_slope_from_dem(dem_path, lon_click, lat_click) if dem_path and Path(dem_path).exists() else None

        st.markdown(f"**Location:** {lat_click:.6f}, {lon_click:.6f}")
        st.markdown(f"- Infiltration: **{infil_val:.3f}**" if infil_val is not None else "- Infiltration: n/a")
        st.markdown(f"- AWC (raw): **{awc_val:.3f}**" if awc_val is not None else "- AWC (raw): n/a")
        st.markdown(f"- Slope (norm ~[0–1]): **{slope_val:.3f}**" if slope_val is not None else "- Slope: n/a")
    else:
        st.write("Click on the map to inspect values at a point.")

st.success("Map ready. Adjust opacity/palette, toggle AOI, and inspect values.")
