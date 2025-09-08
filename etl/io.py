"""
Thin vector I/O adapter using pyogrio when available.
Keeps GeoPandas calls centralized so we can swap engines or defaults later.
"""
from __future__ import annotations
from pathlib import Path
from typing import Optional
import geopandas as gpd

_PYOGRIO_ENGINE = "pyogrio"  # fallbacks to fiona automatically if missing


def read_vector(path: str | Path) -> gpd.GeoDataFrame:
    """Read any OGR vector (GeoPackage, GeoJSON, Shapefile)."""
    return gpd.read_file(path, engine=_PYOGRIO_ENGINE)


def write_vector(
    gdf: gpd.GeoDataFrame,
    path: str | Path,
    driver: Optional[str] = None,
    layer: Optional[str] = None,
    index: bool = False,
) -> None:
    """Write vector with pyogrio engine (defaults sensible for GPKG)."""
    path = Path(path)
    kwargs = {"engine": _PYOGRIO_ENGINE, "index": index}
    if driver is None:
        # Infer driver from extension if possible, else let GeoPandas decide
        ext = path.suffix.lower()
        if ext in {".gpkg", ".geopackage"}:
            driver = "GPKG"
        elif ext in {".geojson", ".json"}:
            driver = "GeoJSON"
    if driver:
        kwargs["driver"] = driver
    if layer:
        kwargs["layer"] = layer
    gdf.to_file(path, **kwargs)


def to_geopackage(
    gdf: gpd.GeoDataFrame,
    gpkg_path: str | Path,
    layer: str,
    index: bool = False,
) -> None:
    """Append/overwrite a layer inside a GeoPackage."""
    write_vector(gdf, gpkg_path, driver="GPKG", layer=layer, index=index)
