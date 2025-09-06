import re, requests, time
from pathlib import Path

try:
    import gdown  # handles Google Drive large-file confirmation pages
except Exception:
    gdown = None

_ID_RE = re.compile(r"/d/([A-Za-z0-9_-]{20,})/|id=([A-Za-z0-9_-]{20,})")

def extract_id(url: str):
    m=_ID_RE.search(url or ""); return (m.group(1) or m.group(2)) if m else None

def is_drive_url(url: str) -> bool:
    return extract_id(url) is not None

def drive_share_to_download(url: str)->str:
    fid=extract_id(url or ""); return f"https://drive.google.com/uc?export=download&id={fid}" if fid else url

def local_target_for(url: str, base_dir: Path, suggested_name: str|None=None)->Path:
    base_dir.mkdir(parents=True, exist_ok=True); fid=extract_id(url or "") or "unknown"; name=suggested_name or f"{fid}.bin"; return base_dir/name

def is_fresh(path: Path, ttl_hours:int=24)->bool:
    if not path.exists(): return False
    return (time.time()-path.stat().st_mtime)/3600.0 <= ttl_hours

def _download_via_requests(path: Path, url: str, chunk=2**20):
    path.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        with open(path, "wb") as f:
            for b in r.iter_content(chunk_size=chunk):
                if b: f.write(b)
    return path

def _download_via_gdown_id(path: Path, file_id: str):
    if gdown is None:
        raise RuntimeError("gdown is not installed; cannot handle Google Drive large-file downloads.")
    path.parent.mkdir(parents=True, exist_ok=True)
    # gdown returns the output path on success, None on failure
    out = gdown.download(id=file_id, output=str(path), quiet=True, fuzzy=True)
    if out is None or not Path(out).exists():
        raise RuntimeError("gdown failed to download the Google Drive file. Check sharing permissions.")
    return Path(out)

def download_to(path: Path, url: str, chunk=2**20):
    """Download URL to path. If it's a Google Drive link, use gdown for robustness."""
    fid = extract_id(url or "")
    if fid:
        return _download_via_gdown_id(path, fid)
    else:
        return _download_via_requests(path, url, chunk=chunk)
