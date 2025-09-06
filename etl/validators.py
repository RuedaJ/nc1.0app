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

_DRIVE_ID_RE = re.compile(r"/d/([A-Za-z0-9_-]{20,})/|id=([A-Za-z0-9_-]{20,})")

def extract_drive_id(url: str) -> str | None:
    m = _DRIVE_ID_RE.search(url or "")
    if not m:
        return None
    return m.group(1) or m.group(2)

def looks_like_drive_url(url: str) -> bool:
    return extract_drive_id(url) is not None
