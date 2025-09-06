import re

def valid_lat(lat_str: str) -> bool:
    try:
        v = float(lat_str)
        return -90.0 <= v <= 90.0
    except Exception:
        return False

def valid_lon(lon_str: str) -> bool:
    try:
        v = float(lon_str)
        return -180.0 <= v <= 180.0
    except Exception:
        return False

# Google Drive file id extractor
_DRIVE_ID_RE = re.compile(r"/d/([A-Za-z0-9_-]{20,})/|id=([A-Za-z0-9_-]{20,})")

def extract_drive_id(url: str) -> str | None:
    m = _DRIVE_ID_RE.search(url or "")
    return (m.group(1) or m.group(2)) if m else None

def looks_like_drive_url(url: str) -> bool:
    return extract_drive_id(url) is not None

def sniff_signature(path: str) -> str:
    """Return a short signature: 'TIFF', 'SQLITE', 'HTML', or 'UNKNOWN'."""
    try:
        with open(path, "rb") as f:
            head = f.read(64)
        if head.startswith(b"II*\x00") or head.startswith(b"MM\x00*"):
            return "TIFF"
        if head.startswith(b"SQLite format 3\x00"):
            return "SQLITE"
        if head.lstrip().lower().startswith(b"<!doctype html") or b"<html" in head.lower():
            return "HTML"
        return "UNKNOWN"
    except Exception:
        return "UNKNOWN"

def validate_geotiff(path: str) -> tuple[bool, str]:
    """Check header looks like TIFF and that rasterio can open it."""
    sig = sniff_signature(path)
    if sig != "TIFF":
        return False, f"Wrong file signature: {sig}"
    try:
        import rasterio as rio
        with rio.open(path) as ds:
            _ = (ds.count, ds.width, ds.height)
        return True, "OK"
    except Exception as e:
        return False, f"{type(e).__name__}: {e}"
