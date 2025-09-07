# -*- coding: utf-8 -*-
# app/pages/01_Mapa.py
import sys
import os
from pathlib import Path

# Ensure repo root on path
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
from pyproj import Transformer

# Map libs
import folium
try:
    from streamlit_folium import st_folium
    HAS_STF = True
except Exception:
    HAS_STF = False

# Optional fast tiling
try:
    from leafmap.common import get_local_tile_layer
    HAS_LTS = True
except Exception:
    HAS_LTS = False

# Fallback colorizing / image embed
import matplotlib.cm as cm
from PIL import Image
import base64
import io

from state import get_state

st.set_page_config(page_title="Mapa - sustai-geo-app", layout="wide")
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
    show_aoi = st.checkbox("Show AOI", value=True)
    show_site = st.checkbox("Show site marker", value=True)
    if not HAS_STF:
        st.info("For click-to-inspect, add streamlit-folium to requirements. Using fallback display.")

# ---------- Helpers ----------
def add_aoi_to_folium_map(fm, aoi_path):
    try:
        gdf = gpd.read_file(aoi_path)
        gj = folium.GeoJson(
            gdf.to_crs(4326).__geo_interface__,
            name="AOI",
            style_function=lambda _:
                {"color": "#0066CC", "weight": 2, "fill": False},
        )
        gj.add_to(fm)
        minx, miny, maxx, maxy = gdf.to_crs(4326).total_bounds
        fm.fit_bounds([[miny, minx], [maxy, maxx]])
    except Exception as e:
        st.warning("Could not render AOI: {}".format(e))

def add_site_marker_to_folium_map(fm, lat, lon):
    folium.CircleMarker(
        location=[lat, lon],
        radius=5,
        color="#DC3545",
        fill=True,
        fill_opacity=1,
        popup="Site",
    ).add_to(fm)

def _rgba_data_uri(rgba):
    im = Image.fromarray(rgba, mode="RGBA")
    buf = io.BytesIO()
    im.save(buf, format="PNG")
    b64 = base64.b64encode(buf.getvalue()).decode("ascii")
    return "data:image/png;base64,{}".format(b64)

@st.cache_data(show_spinner=False)
def _prep_overlay_for_display(tif_path, max_px=1024):
    """
    Load -> reproject to EPSG:4326 -> optional downsample.
    Returns (arr_float32, (south, west, north, east), nodata_value).
    """
    with rio.open(tif_path) as src:
        nodata = src.nodata
        st.caption(
            "Infiltration GeoTIFF — shape={}x{}, crs={}, dtype={}, nodata={}".format(
                src.height, src.width, src.crs, src.dtypes[0], nodata
            )
        )
    da = rxr.open_rasterio(tif_path, masked=True).squeeze()
    if not da.rio.crs:
        raise ValueError("Infiltration raster has no CRS.")
    da_4326 = da.rio.reproject("EPSG:4326", resampling=Resampling.bilinear)

    # Downsample to keep frontend fast
    max_px = int(max(256, max_px))
    h = da_4326.sizes[da_4326.dims[0]]
    w = da_4326.sizes[da_4326.dims[1]]
    if max(h, w) > max_px:
        resx, resy = da_4326.rio.resolution()
        scale = float(max_px) / float(max(h, w))
        da_4326 = da_4326.rio.reproject(
            "EPSG:4326",
            resolution=(abs(resx) / scale, abs(resy) / scale),
            resampling=Resampling.bilinear,
        )

    arr = da_4326.values.astype("float32")
    bounds = da_4326.rio.bounds()
    return arr, bounds, nodata

def add_raster_to_folium(fm, tif_path, name, opacity):
    """
    Try fast tiles (localtileserver); else cached ImageOverlay with:
    - nodata-aware transparency
    - percentile stretch (2–98%) for good contrast
    """
    if HAS_LTS:
        try:
            tile_layer, _client = get_local_tile_layer(
                tif_path,
                name=name,
                palette="Blues",
                vmin=0, vmax=1,
                opacity=opacity,
            )
            tile_layer.add_to(fm)
            folium.LayerControl(collapsed=False).add_to(fm)
            return
        except Exception as e:
            st.info("localtileserver path unavailable ({}). Using ImageOverlay fallback.".format(e))

    # Fallback overlay
    try:
        arr, (south, west, north, east), nodata = _prep_overlay_for_display(tif_path)
    except Exception as e:
        st.error("Failed to prepare infiltration overlay: {}".format(e))
        return

    if arr.size == 0 or arr.shape[0] == 0 or arr.shape[1] == 0:
        st.warning("Infiltration raster has zero-size after clipping/reprojection. Check AOI overlap and input coverage.")
        return

    finite = np.isfinite(arr)
    valid = finite.copy()
    if nodata is not None and np.isfinite(nodata):
        valid &= arr != float(nodata)

    valid_ratio = float(valid.sum()) / float(valid.size) if valid.size else 0.0
    st.caption("Infiltration display: valid%={:.1f}".format(100.0 * valid_ratio))

    if valid.any():
        # Percentile stretch for robust contrast
        qlo = float(np.nanpercentile(arr[valid], 2))
        qhi = float(np.nanpercentile(arr[valid], 98))
        if qhi > qlo:
            arr_norm = (arr - qlo) / (qhi - qlo)
        else:
            arr_norm = np.zeros_like(arr, dtype="float32")
    else:
        arr_norm = np.zeros_like(arr, dtype="float32")

    rgba = (cm.get_cmap("Blues")(np.clip(arr_norm, 0, 1)) * 255).astype("uint8")
    # Make invalid pixels transparent and apply global opacity
    rgba[~valid, 3] = 0
    rgba[..., 3] = (rgba[..., 3].astype("float32") * float(opacity)).astype("uint8")

    folium.raster_layers.ImageOverlay(
        image=_rgba_data_uri(rgba),
        bounds=[[south, west], [north, east]],
        name=name,
        opacity=1.0,  # alpha baked in
        interactive=False,
        cross_origin=False,
        zindex=500,
    ).add_to(fm)
    folium.LayerControl(collapsed=False).add_to(fm)

def sample_tiff_value(tif_path, lon, lat):
    try:
        with rio.open(tif_path) as ds:
            if ds.crs is None:
                return None
            transformer = Transformer.from_crs("EPSG:4326", ds.crs, always_xy=True)
            x, y = transformer.transform(lon, lat)
            val = list(ds.sample([(x, y)]))[0][0]
            if isinstance(val, np.ndarray):
                val = float(val[0])
            if val is None:
                return None
            val = float(val)
            if np.isnan(val):
                return None
            return val
    except Exception:
        return None

def sample_slope_from_dem(dem_path, lon, lat):
    """Approximate normalized slope [0,1] from a 3x3 window at click location."""
    try:
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
            a, b, c_aff, d, e, f_aff, _, _, _ = ds.transform
            gy, gx = np.gradient(arr)
            with np.errstate(invalid="ignore", divide="ignore"):
                slope_rad = np.arctan(np.sqrt((gx * (1.0 / abs(a))) ** 2 + (gy * (1.0 / abs(e))) ** 2))
            slope_deg = np.degrees(slope_rad)
            if np.isnan(slope_deg).all():
                return None
            smin = float(np.nanmin(slope_deg))
            smax = float(np.nanmax(slope_deg))
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

# AOI and site
if aoi_geojson and Path(aoi_geojson).exists() and show_aoi:
    add_aoi_to_folium_map(fm, aoi_geojson)
if show_site and coords:
    add_site_marker_to_folium_map(fm, coords[0], coords[1])

# Raster layer (tiles if possible; otherwise embedded PNG overlay)
add_raster_to_folium(fm, str(infil_path), "Infiltration", opacity)

# Simple legend
legend_html = """
<div style="position:relative;margin-top:6px;">
  <div style="font:12px/1.2 sans-serif;color:#333;">Infiltration (0 -> 1)</div>
  <div style="height:10px;background:linear-gradient(to right,#f7fbff,#deebf7,#9ecae1,#4292c6,#08519c);"></div>
  <div style="display:flex;justify-content:space-between;font:11px sans-serif;color:#555;">
    <span>0.0</span><span>0.5</span><span>1.0</span>
  </div>
</div>
"""
st.markdown(legend_html, unsafe_allow_html=True)

if HAS_STF:
    st_data = st_folium(fm, height=640, returned_objects=["last_clicked"])
    last = st_data.get("last_clicked")
    with st.expander("Inspector", expanded=bool(last)):
        if last:
            lat_click = last["lat"]
            lon_click = last["lng"]
            infil_val = sample_tiff_value(infil_path, lon_click, lat_click)
            awc_path = S.get("awc_path")
            awc_val = sample_tiff_value(awc_path, lon_click, lat_click) if awc_path and Path(awc_path).exists() else None
            dem_path = S.get("dem_path")
            slope_val = sample_slope_from_dem(dem_path, lon_click, lat_click) if dem_path and Path(dem_path).exists() else None

            st.markdown("**Location:** {:.6f}, {:.6f}".format(lat_click, lon_click))
            st.markdown("- Infiltration: **{:.3f}**".format(infil_val) if infil_val is not None else "- Infiltration: n/a")
            st.markdown("- AWC (raw): **{:.3f}**".format(awc_val) if awc_val is not None else "- AWC (raw): n/a")
            st.markdown("- Slope (norm ~[0-1]): **{:.3f}**".format(slope_val) if slope_val is not None else "- Slope: n/a")
        else:
            st.write("Click on the map to inspect values at a point.")
else:
    # Static snapshot + manual sampler
    st.info("Install streamlit-folium for click inspector. Showing map and manual sampler.")
    html_path = Path("data/processed/_tmp_map.html")
    html_path.parent.mkdir(parents=True, exist_ok=True)
    fm.save(str(html_path))
    st.components.v1.iframe(src=str(html_path), height=640)

    with st.expander("Inspector (fallback: manual point)"):
        col1, col2, _ = st.columns([1, 1, 2])
        with col1:
            lat_m = st.number_input("Lat", value=float(coords[0]) if coords else 41.65, format="%.6f")
        with col2:
            lon_m = st.number_input("Lon", value=float(coords[1]) if coords else -0.89, format="%.6f")
        if st.button("Sample at point"):
            infil_val = sample_tiff_value(infil_path, lon_m, lat_m)
            awc_path = S.get("awc_path")
            dem_path = S.get("dem_path")
            awc_val = sample_tiff_value(awc_path, lon_m, lat_m) if awc_path and Path(awc_path).exists() else None
            slope_val = sample_slope_from_dem(dem_path, lon_m, lat_m) if dem_path and Path(dem_path).exists() else None
            st.markdown("**Location:** {:.6f}, {:.6f}".format(lat_m, lon_m))
            st.markdown("- Infiltration: **{:.3f}**".format(infil_val) if infil_val is not None else "- Infiltration: n/a")
            st.markdown("- AWC (raw): **{:.3f}**".format(awc_val) if awc_val is not None else "- AWC (raw): n/a")
            st.markdown("- Slope (norm ~[0-1]): **{:.3f}**".format(slope_val) if slope_val is not None else "- Slope: n/a")

st.success("Map ready. Adjust opacity, toggle AOI, and inspect values.")
