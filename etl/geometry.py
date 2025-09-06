from pathlib import Path
import geopandas as gpd
from shapely.geometry import Point

def _utm_for(lon: float, lat: float) -> str:
    zone = int((lon + 180) // 6) + 1
    epsg = 32600 + zone if lat >= 0 else 32700 + zone
    return f"EPSG:{epsg}"

def aoi_from_latlon(lat: float, lon: float, radius_km: float, out_crs: str) -> gpd.GeoDataFrame:
    """Build a circular AOI buffered in meters using a local UTM CRS, then return in out_crs."""
    gdf = gpd.GeoDataFrame(geometry=[Point(lon, lat)], crs="EPSG:4326")
    utm = _utm_for(lon, lat)
    gdf_utm = gdf.to_crs(utm)
    aoi_utm = gdf_utm.buffer(radius_km * 1000.0)
    aoi = gpd.GeoDataFrame(geometry=aoi_utm, crs=utm).to_crs(out_crs)
    return aoi

def save_geojson(gdf: gpd.GeoDataFrame, path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_file(path, driver="GeoJSON")
