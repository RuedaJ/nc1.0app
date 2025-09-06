import re
def valid_lat(lat_str: str) -> bool:
    try: v=float(lat_str); return -90<=v<=90
    except: return False
def valid_lon(lon_str: str) -> bool:
    try: v=float(lon_str); return -180<=v<=180
    except: return False
_DRIVE_ID_RE = re.compile(r"/d/([A-Za-z0-9_-]{20,})/|id=([A-Za-z0-9_-]{20,})")
def extract_drive_id(url: str) -> str | None:
    m = _DRIVE_ID_RE.search(url or ""); 
    return (m.group(1) or m.group(2)) if m else None
def looks_like_drive_url(url: str) -> bool:
    return extract_drive_id(url) is not None

def validate_geotiff(path: str) -> tuple[bool, str]:
    """Try to open with rasterio; return (ok, message)."""
    try:
        import rasterio as rio
        with rio.open(path) as ds:
            _ = ds.count; _ = ds.width; _ = ds.height
        return True, "OK"
    except Exception as e:
        return False, f"{type(e).__name__}: {e}"
