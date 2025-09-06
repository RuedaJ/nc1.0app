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
    """Return a short signature: 'TIFF', 'BIGTIFF', 'VRT', 'SQLITE', 'HTML', 'ZIP', 'GZIP', 'PDF', 'JSON', 'UNKNOWN'."""
    try:
        with open(path, "rb") as f:
            head = f.read(1024)
        h = head.lstrip()
        # TIFF classic
        if head.startswith(b"II*\x00") or head.startswith(b"MM\x00*"):
            return "TIFF"
        # BigTIFF (0x2B)
        if head.startswith(b"II+\x00") or head.startswith(b"MM\x00+"):
            return "BIGTIFF"
        # ZIP (often zipped GeoTIFF)
        if head.startswith(b"PK\x03\x04") or head.startswith(b"PK\x05\x06") or head.startswith(b"PK\x07\x08"):
            return "ZIP"
        # GZIP
        if head.startswith(b"\x1f\x8b"):
            return "GZIP"
        # PDF
        if head.startswith(b"%PDF-"):
            return "PDF"
        # SQLite (GeoPackage)
        if head.startswith(b"SQLite format 3\x00"):
            return "SQLITE"
        # VRT (XML)
        if h.startswith(b"<VRTDataset") or h.startswith(b"<?xml"):
            return "VRT"
        # HTML (loose)
        low = head.lower()
        if b"<html" in low or h.startswith(b"<!doctype html"):
            return "HTML"
        # JSON
        if h.startswith(b"{") or h.startswith(b"["):
            return "JSON"
        return "UNKNOWN"
    except Exception:
        return "UNKNOWN"

def validate_geotiff(path: str) -> tuple[bool, str]:
    sig = sniff_signature(path)
    if sig in ("TIFF", "BIGTIFF", "VRT"):
        try:
            import rasterio as rio
            with rio.open(path) as ds:
                _ = (ds.count, ds.width, ds.height)
            return True, "OK" if sig != "VRT" else "OK (VRT)"
        except Exception as e:
            return False, f"{type(e).__name__}: {e}"
    return False, f"Wrong file signature: {sig}"
