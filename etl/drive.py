import os, re, requests, time, shutil
from pathlib import Path
from typing import Optional, Tuple

try:
    import gdown  # robust Google Drive downloader
except Exception:
    gdown = None

# --- Patterns ---
FILE_RE   = re.compile(r"/file/d/([A-Za-z0-9_-]{20,})")
ID_RE     = re.compile(r"[?&]id=([A-Za-z0-9_-]{20,})")
FOLDER_RE = re.compile(r"/folders/([A-Za-z0-9_-]{20,})")

RASTER_EXTS = (".tif", ".tiff", ".vrt", ".img")  # extend as needed

# ---------------- Core parsing ----------------
def parse_drive_url(url: str) -> Tuple[Optional[str], Optional[str]]:
    """Return (kind, id) where kind in {'file','folder'} or (None,None)."""
    if not url:
        return None, None
    m = FOLDER_RE.search(url)
    if m:
        return "folder", m.group(1)
    m = FILE_RE.search(url)
    if m:
        return "file", m.group(1)
    m = ID_RE.search(url)
    if m:
        # assume file when id= is provided
        return "file", m.group(1)
    return None, None

def is_drive_url(url: str) -> bool:
    k, i = parse_drive_url(url or "")
    return k is not None and i is not None

# ---------------- Local caching helpers ----------------
def local_target_for(url: str, base_dir: Path, suggested_name: str | None = None) -> Path:
    base_dir.mkdir(parents=True, exist_ok=True)
    kind, fid = parse_drive_url(url or "")
    name = suggested_name or (f"{fid}.bin" if fid else "download.bin")
    return base_dir / name

def is_fresh(path: Path, ttl_hours: int = 24) -> bool:
    if not path.exists():
        return False
    return (time.time() - path.stat().st_mtime) / 3600.0 <= ttl_hours

def _download_via_requests(path: Path, url: str, chunk=2**20) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        with open(path, "wb") as f:
            for b in r.iter_content(chunk_size=chunk):
                if b:
                    f.write(b)
    return path

def download_to(path: Path, url: str, chunk=2**20) -> Path:
    if not url:
        raise ValueError("Empty URL for download.")
    path.parent.mkdir(parents=True, exist_ok=True)
    return _download_via_requests(path, url, chunk=chunk)

# ---------------- Folder-aware Drive download ----------------
def _find_first_raster_in_dir(d: Path, exts=RASTER_EXTS) -> Path | None:
    if not d.exists() or not d.is_dir():
        return None
    for ext in exts:
        cands = sorted(d.rglob(f"*{ext}"))  # search recursively for nested files
        if cands:
            return cands[0]
    return None

def _cleanup_parts(d: Path) -> None:
    if not d.exists() or not d.is_dir():
        return
    for p in d.glob("*.part"):
        try:
            p.unlink()
        except Exception:
            pass

def download_drive_any(url: str, dest_dir: Path, prefer_exts=RASTER_EXTS) -> Path:
    """Download a Google Drive file or folder; return a concrete raster file path.
    - For *folder* links: downloads folder into dest_dir, returns first raster file.
    - For *file* links: downloads into dest_dir, returns the file path (or the first raster found).
    """
    if gdown is None:
        raise RuntimeError("gdown is not installed; cannot download Google Drive files.")

    kind, fid = parse_drive_url(url or "")
    if not fid:
        raise ValueError("Not a valid Google Drive link: missing id.")

    dest_dir.mkdir(parents=True, exist_ok=True)
    _cleanup_parts(dest_dir)

    if kind == "folder":
        out = gdown.download_folder(id=fid, output=str(dest_dir), quiet=True, use_cookies=False)
        if not out:
            raise RuntimeError("gdown.download_folder returned no files. Check sharing permissions.")
        # try to find a raster inside the downloaded folder tree
        ras = _find_first_raster_in_dir(dest_dir, prefer_exts)
        if ras is None:
            raise RuntimeError("Downloaded folder contains no .tif/.tiff/.vrt/.img. Place a raster in the folder or link directly to the file.")
        return ras

    # kind == "file"
    out = gdown.download(id=fid, output=str(dest_dir), quiet=True, fuzzy=True, use_cookies=False)
    if out is None:
        raise RuntimeError("gdown failed to download the Drive file. Check that sharing is 'Anyone with the link'.")
    outp = Path(out)
    if outp.is_dir():
        cand = _find_first_raster_in_dir(outp, prefer_exts) or _find_first_raster_in_dir(dest_dir, prefer_exts)
        if cand:
            return cand
    if outp.suffix.lower() in prefer_exts and outp.exists():
        return outp

    cand = _find_first_raster_in_dir(dest_dir, prefer_exts)
    if cand:
        return cand

    raise RuntimeError(f"Downloaded file is not a raster: {outp.name}. Provide a direct .tif link or a Drive *folder* containing the .tif.")

# ---------------- Backward-compat shims ----------------
def extract_id(url: str) -> str | None:
    """Legacy: return Drive file/folder id (prefers whichever is present)."""
    _, fid = parse_drive_url(url or "")
    return fid

def download_drive_file(url: str, dest_dir: Path) -> Path:
    """Legacy: behave like the old function by delegating to folder-aware downloader."""
    return download_drive_any(url, dest_dir)
