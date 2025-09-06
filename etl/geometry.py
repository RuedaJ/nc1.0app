from pathlib import Path
import geopandas as gpd
from shapely.geometry import Point

WGS84 = "EPSG:4326"

def aoi_from_latlon(lat: float, lon: float, radius_km: float, crs_out: str = "EPSG:25830") -> gpd.GeoDataFrame:
    g = gpd.GeoSeries([Point(float(lon), float(lat))], crs=WGS84)
    g_out = g.to_crs(crs_out)
    aoi = g_out.buffer(radius_km * 1000.0)
    return gpd.GeoDataFrame(geometry=aoi, crs=crs_out)

def save_geojson(gdf: gpd.GeoDataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_file(path, driver="GeoJSON")
    return path
