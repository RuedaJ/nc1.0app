"""
Local tile server helper for raster layers (COG/GeoTIFF).
"""
from __future__ import annotations
from pathlib import Path

try:
    from localtileserver import TileClient, get_leaflet_tile_layer
except Exception as e:  # pragma: no cover
    TileClient = None
    get_leaflet_tile_layer = None
    _IMPORT_ERROR = e
else:
    _IMPORT_ERROR = None


def raster_tile_layer(path: str | Path):
    """
    Return a Folium-compatible tile layer and the TileClient for a raster file.
    Raises if localtileserver isn't installed.
    """
    if TileClient is None:
        raise RuntimeError(
            "localtileserver is not available. Install it or disable tiled rasters."
        ) from _IMPORT_ERROR
    client = TileClient(str(Path(path)))
    tile = get_leaflet_tile_layer(
        client,
        attribution="© data sources",
        name=Path(path).name,
        opacity=1.0,
        shown=True,
    )
    return tile, client
