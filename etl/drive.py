import re, requests, time
from pathlib import Path

_ID_RE = re.compile(r"/d/([A-Za-z0-9_-]{20,})/|id=([A-Za-z0-9_-]{20,})")

def extract_id(url: str) -> str | None:
    m = _ID_RE.search(url or "")
    if not m:
        return None
    return m.group(1) or m.group(2)

def drive_share_to_download(url: str) -> str:
    file_id = extract_id(url or "")
    if file_id:
        return f"https://drive.google.com/uc?export=download&id={file_id}"
    return url

def local_target_for(url: str, base_dir: Path, suggested_name: str | None = None) -> Path:
    base_dir.mkdir(parents=True, exist_ok=True)
    file_id = extract_id(url or "") or "unknown"
    name = suggested_name or f"{file_id}.bin"
    return base_dir / name

def is_fresh(path: Path, ttl_hours: int = 24) -> bool:
    if not path.exists():
        return False
    import time
    age_hours = (time.time() - path.stat().st_mtime) / 3600.0
    return age_hours <= ttl_hours

def download_to(path: Path, url: str, chunk=2**20):
    path.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        with open(path, "wb") as f:
            for b in r.iter_content(chunk_size=chunk):
                if b:
                    f.write(b)
    return path
