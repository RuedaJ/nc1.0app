# app/pages/01_Mapa.py
import sys, os
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
from shapely.geometry import box
from pyproj import Transformer

# Map libs
import folium
try:
    from streamlit_folium import st_folium
    HAS_STF = True
except Exception:
    HAS_STF = False

# We still use leafmap for its helpers/visual convenience on the fallback path
import leafmap.foliumap as leafmap
import matplotlib.cm as cm  # for ImageOverlay colorizing

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
    show_aoi = st.checkbox("Show AOI", value=True)
    show_site = st.checkbox("Show site marker", value=True)
    if not HAS_STF:
        st.info("For click-to-inspect, add `streamlit-folium` to requirements. Using fallback display.")

# ---------- Shared helpers ----------
def add_aoi_to_folium_map(fm: folium.Map, aoi_path: str):
    try:
        gdf = gpd.read_file(aoi_path)
        gj = folium.GeoJson(gdf.to_crs(4326).__geo_interface__,
                            name="AOI",
                            style_function=lambda _:
                                {"color": "#0066CC", "weight": 2, "fill": False})
        gj.add_to(fm)
        # Fit to AOI bounds
        minx, miny, maxx, maxy = gdf.to_crs(4326).total_bounds
        fm.fit_bounds([[miny, minx], [maxy, maxx]])
    except Exception as e:
        st.warning(f"Could not render AOI: {e}")

def add_site_marker_to_folium_map(fm: folium.Map, lat: float, lon: float):
    folium.CircleMarker(location=[lat, lon], radius=5,
                        color="#DC3545", fill=True, fill_opacity=1,
                        popup="Site").add_to(fm)

def add_raster_tiles_to_folium(fm: folium.Map, tif_path: str, name: str, opacity: float):
    """Use localtileserver via leafmap.common if available; else raise ImportError."""
    from leafmap.common import get_local_tile_layer  # may raise ImportError
    tile_layer, _client = get_local_tile_layer(
        tif_path,
        name=name,
        palette="Blues",
        vmin=0, vmax=1,
        opacity=opacity,
    )
    tile_layer.add_to(fm)
    folium.LayerControl(collapsed=False).add_to(fm)

def add_raster_imageoverlay_to_folium(fm: folium.Map, tif_path: str, name: str, opacity: float):
    """Serverless fallback: reproject to EPSG:4326, downsample, colorize with matplotlib, add as ImageOverlay."""
    da = rxr.open_rasterio(tif_path, masked=True).squeeze()
    da_4326 = da.rio.reproject("EPSG:4326", resampling=Resampling.bilinear)

    # Downsample if large
    max_px = 1024
    h, w = da_4326.sizes[da_4326.dims[0]], da_4326.sizes[da_4326.dims[1]]
    scale = min(1.0, max_px / max(h, w))
    if scale < 1.0:
        resx, resy = da_4326.rio.resolution()
        da_4326 = da_4326.rio.reproject(
            "EPSG:4326",
            resolution=(abs(resx) / scale, abs(resy) / scale),
            resampling=Resampling.bilinear,
        )

    arr = da_4326.values.astype("float32")
    finite = np.isfinite(arr)
    if finite.any():
        vmin, vmax = np.nanmin(arr[finite]), np.nanmax(arr[finite])
        if vmax > vmin:
            arr = (arr - vmin) / (vmax - vmin)
        else:
            arr = np.zeros_like(arr)
    else:
        arr = np.zeros_like(arr)

    rgba = (cm.get_cmap("Blues")(np.clip(arr, 0, 1)) * 255).astype("uint8")
    south, west, north, east = da_4326.rio.bounds()
    folium.raster_layers.ImageOverlay(
        image=rgba,
        bounds=[[south, west], [north, east]],
        name=name,
        opacity=opacity,
        interactive=False,
        cross_origin=False,
        zindex=1,
    ).add_to(fm)
    folium.LayerControl(collapsed=False).add_to(fm)

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
        with rio.open(dem_path) as ds:
            if ds.crs is None:
                return None
            transformer = Transformer.from_crs("EPSG:4326", ds.crs, always_xy=True)
            x, y = transformer.transform(lon, lat)
            r, c = ds.index(x, y)
            win = Window(col_off=max(c-1,0), row_off=max(r-1,0),
                         width=3 if c+2 < ds.width  else min(3, ds.width - max(c-1,0)),
                         height=3 if r+2 < ds.height else min(3, ds.height - max(r-1,0)))
            arr = ds.read(1, window=win, boundless=True, fill_value=np.nan).astype("float32")
            if np.isnan(arr).all():
                return None
            a, b, c_aff, d, e, f_aff, _, _, _ = ds.transform
            gy, gx = np.gradient(arr)
            with np.errstate(invalid="ignore", divide="ignore"):
                slope_rad = np.arctan(np.sqrt((gx * (1.0/abs(a)))**2 + (gy * (1.0/abs(e)))**2))
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
if HAS_STF:
    # Build a pure folium map so click events work
    base_loc = (coords[0], coords[1]) if coords else (41.65, -0.89)
    fm = folium.Map(location=base_loc, zoom_start=12, control_scale=True, tiles="CartoDB positron")

    # AOI and site
    if aoi_geojson and Path(aoi_geojson).exists() and show_aoi:
        add_aoi_to_folium_map(fm, aoi_geojson)
    if show_site and coords:
        add_site_marker_to_folium_map(fm, coords[0], coords[1])

    # Raster layer
    try:
        add_raster_tiles_to_folium(fm, str(infil_path), "Infiltration", opacity)
    except Exception as e:
        st.info(f"`localtileserver` path unavailable ({e}). Using ImageOverlay fallback.")
        add_raster_imageoverlay_to_folium(fm, str(infil_path), "Infiltration", opacity)

    # Show + click inspector
    st_data = st_folium(fm, height=640, returned_objects=["last_clicked"])
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
else:
    # Fallback display without streamlit_folium: use leafmap and manual sampler
    base_loc = (coords[0], coords[1]) if coords else (41.65, -0.89)
    m = leafmap.Map(center=base_loc, zoom=12, draw_control=False, measure_control=False)

    # AOI and site
    if aoi_geojson and Path(aoi_geojson).exists() and show_aoi:
        try:
            gdf = gpd.read_file(aoi_geojson)
            m.add_gdf(gdf, layer_name="AOI", style={"color": "#0066CC", "weight": 2, "fill": False})
            minx, miny, maxx, maxy = gdf.to_crs(4326).total_bounds
            m.fit_bounds([[miny, minx], [maxy, maxx]])
        except Exception as e:
            st.warning(f"Could not render AOI: {e}")
    if show_site and coords:
        # leafmap convenience
        m.add_marker(location=(coords[0], coords[1]), popup="Site")

    # Raster layer (ImageOverlay fallback)
    try:
        # Try tiling via localtileserver on leafmap (if installed)
        m.add_raster(str(infil_path), palette="Blues", layer_name="Infiltration", opacity=opacity, vmin=0, vmax=1)
    except Exception:
        # Fallback: ImageOverlay style (use folium layer added via add_layer)
        da = rxr.open_rasterio(str(infil_path), masked=True).squeeze()
        da_4326 = da.rio.reproject("EPSG:4326", resampling=Resampling.bilinear)
        max_px = 1024
        h, w = da_4326.sizes[da_4326.dims[0]], da_4326.sizes[da_4326.dims[1]]
        scale = min(1.0, max_px / max(h, w))
        if scale < 1.0:
            resx, resy = da_4326.rio.resolution()
            da_4326 = da_4326.rio.reproject(
                "EPSG:4326",
                resolution=(abs(resx) / scale, abs(resy) / scale),
                resampling=Resampling.bilinear,
            )
        arr = da_4326.values.astype("float32")
        finite = np.isfinite(arr)
        if finite.any():
            vmin, vmax = np.nanmin(arr[finite]), np.nanmax(arr[finite])
            if vmax > vmin:
                arr = (arr - vmin) / (vmax - vmin)
            else:
                arr = np.zeros_like(arr)
        else:
            arr = np.zeros_like(arr)
        rgba = (cm.get_cmap("Blues")(np.clip(arr, 0, 1)) * 255).astype("uint8")
        south, west, north, east = da_4326.rio.bounds()
        overlay = folium.raster_layers.ImageOverlay(
            image=rgba,
            bounds=[[south, west], [north, east]],
            name="Infiltration",
            opacity=opacity,
            interactive=False,
            cross_origin=False,
            zindex=1,
        )
        m.add_layer(overlay)

    m.to_streamlit(height=640)

    # Manual inspector
    with st.expander("Inspector (fallback: manual point)"):
        col1, col2, _ = st.columns([1,1,2])
        with col1:
            lat_m = st.number_input("Lat", value=float(coords[0]) if coords else 41.65, format="%.6f")
        with col2:
            lon_m = st.number_input("Lon", value=float(coords[1]) if coords else -0.89, format="%.6f")
        if st.button("Sample at point"):
            infil_val = sample_tiff_value(infil_path, lon_m, lat_m)
            awc_path = S.get("awc_path"); dem_path = S.get("dem_path")
            awc_val = sample_tiff_value(awc_path, lon_m, lat_m) if awc_path and Path(awc_path).exists() else None
            slope_val = sample_slope_from_dem(dem_path, lon_m, lat_m) if dem_path and Path(dem_path).exists() else None
            st.markdown(f"**Location:** {lat_m:.6f}, {lon_m:.6f}")
            st.markdown(f"- Infiltration: **{infil_val:.3f}**" if infil_val is not None else "- Infiltration: n/a")
            st.markdown(f"- AWC (raw): **{awc_val:.3f}**" if awc_val is not None else "- AWC (raw): n/a")
            st.markdown(f"- Slope (norm ~[0–1]): **{slope_val:.3f}**" if slope_val is not None else "- Slope: n/a")

st.success("Map ready. Toggle layers and inspect values.")
