import os, re, requests, time, shutil, zipfile
from pathlib import Path
from typing import Optional, Tuple, Iterable

try:
    import gdown  # robust Google Drive downloader
except Exception:
    gdown = None

# ---------------- Signatures ----------------
def _sniff_signature(path: str) -> str:
    try:
        with open(path, "rb") as f:
            head = f.read(2048)
        h = head.lstrip()
        if head.startswith(b"II*\x00") or head.startswith(b"MM\x00*"): return "TIFF"
        if head.startswith(b"II+\x00") or head.startswith(b"MM\x00+"): return "BIGTIFF"
        if head.startswith(b"PK\x03\x04") or head.startswith(b"PK\x05\x06") or head.startswith(b"PK\x07\x08"): return "ZIP"
        if head.startswith(b"\x1f\x8b"): return "GZIP"
        if head.startswith(b"%PDF-"): return "PDF"
        if head.startswith(b"SQLite format 3\x00"): return "SQLITE"
        if h.startswith(b"<VRTDataset") or h.startswith(b"<?xml"): return "VRT"
        low = head.lower()
        if b"<html" in low or h.startswith(b"<!doctype html"): return "HTML"
        if h.startswith(b"{") or h.startswith(b"["): return "JSON"
        return "UNKNOWN"
    except Exception:
        return "UNKNOWN"

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

# ---------------- Helpers ----------------
def _cleanup_parts(d: Path) -> None:
    if not d.exists() or not d.is_dir():
        return
    for p in d.glob("*.part"):
        try:
            p.unlink()
        except Exception:
            pass

def _iter_files_recursive(root: Path) -> Iterable[Path]:
    for p in root.rglob("*"):
        if p.is_file():
            yield p

def _find_first_raster_by_ext(d: Path, exts=RASTER_EXTS) -> Path | None:
    if not d.exists() or not d.is_dir():
        return None
    for ext in exts:
        cands = sorted(d.rglob(f"*{ext}"))  # recursive search
        if cands:
            return cands[0]
    return None

def _maybe_unzip(path: Path, into: Path) -> Path | None:
    try:
        if zipfile.is_zipfile(path):
            into.mkdir(parents=True, exist_ok=True)
            with zipfile.ZipFile(path, "r") as z:
                z.extractall(into)
            cand = _find_first_raster_by_ext(into, RASTER_EXTS)
            if cand:
                return cand
            # if none by ext, scan by signature
            for f in _iter_files_recursive(into):
                sig = _sniff_signature(str(f))
                if sig in ("TIFF","BIGTIFF","VRT"):
                    return f
    except Exception:
        return None
    return None

def _scan_dir_by_signature(d: Path) -> Path | None:
    """Recursively scan all files, detect TIFF/BIGTIFF/VRT by signature.
    Unzip ZIPs on the fly into the same directory and continue scanning."""
    for f in _iter_files_recursive(d):
        sig = _sniff_signature(str(f))
        if sig in ("TIFF","BIGTIFF","VRT"):
            return f
        if sig == "ZIP":
            cand = _maybe_unzip(f, d)
            if cand:
                return cand
    return None

# ---------------- Folder-aware Drive download ----------------
def download_drive_any(url: str, dest_dir: Path, prefer_exts=RASTER_EXTS) -> Path:
    """Download a Google Drive file or folder; return a concrete raster file path.
    Handles:
      - Drive *folders* (download_folder) → pick raster by ext or signature.
      - Drive *files* (download) → accept no-extension files via signature.
      - Auto-unzip ZIPs and return first raster.
    """
    if gdown is None:
        raise RuntimeError("gdown is not installed; cannot download Google Drive files.")

    kind, fid = parse_drive_url(url or "")
    if not fid:
        raise ValueError("Not a valid Google Drive link: missing id.")

    dest_dir.mkdir(parents=True, exist_ok=True)
    _cleanup_parts(dest_dir)

    if kind == "folder":
        out_list = gdown.download_folder(id=fid, output=str(dest_dir), quiet=True, use_cookies=True)
        if not out_list:
            raise RuntimeError("gdown.download_folder returned no files. Check sharing permissions.")
        # First try by extension
        ras = _find_first_raster_by_ext(dest_dir, prefer_exts)
        if ras:
            return ras
        # Then try by signature (covers extension-less)
        ras = _scan_dir_by_signature(dest_dir)
        if ras:
            return ras
        raise RuntimeError("Downloaded folder contains no GeoTIFF/VRT. Put a .tif in the folder or link directly to the file.")

    # kind == "file"
    out = gdown.download(id=fid, output=str(dest_dir), quiet=True, fuzzy=True, use_cookies=True)
    if out is None:
        raise RuntimeError("gdown failed to download the Drive file. Check that sharing is 'Anyone with the link'.")
    outp = Path(out)

    if outp.is_dir():
        # Sometimes Drive returns a directory (shortcut/package). Scan it.
        cand = _find_first_raster_by_ext(outp, prefer_exts) or _find_first_raster_by_ext(dest_dir, prefer_exts)
        if cand:
            return cand
        # No ext match; scan by signature
        cand = _scan_dir_by_signature(outp) or _scan_dir_by_signature(dest_dir)
        if cand:
            return cand
        raise RuntimeError("Downloaded item is a directory without rasters. Put a .tif inside or link directly to the .tif.")

    # If known raster extension
    if outp.suffix.lower() in prefer_exts and outp.exists():
        return outp

    # No extension or unknown ext: inspect signature
    sig = _sniff_signature(str(outp))
    if sig in ("TIFF", "BIGTIFF", "VRT"):
        return outp
    if sig == "ZIP":
        cand = _maybe_unzip(outp, dest_dir)
        if cand:
            return cand
        raise RuntimeError("Downloaded ZIP does not contain a recognized raster (.tif/.tiff/.vrt/.img).")
    if sig in ("HTML", "JSON", "PDF", "UNKNOWN"):
        # As a last attempt, scan dest_dir recursively (some Drive cases rename files oddly)
        cand = _scan_dir_by_signature(dest_dir)
        if cand:
            return cand
        raise RuntimeError(f"Downloaded file is not a raster (sig={sig}). Ensure the Drive link points to the actual .tif and sharing is 'Anyone with the link'.")

    # Last resort
    cand = _find_first_raster_by_ext(dest_dir, prefer_exts)
    if cand:
        return cand

    raise RuntimeError(f"Downloaded file is not a raster: {outp.name}. Provide a direct .tif link or a Drive *folder* containing the .tif.")

# ---------------- Backward-compat shims ----------------
def extract_id(url: str) -> str | None:
    _, fid = parse_drive_url(url or "")
    return fid

def download_drive_file(url: str, dest_dir: Path) -> Path:
    return download_drive_any(url, dest_dir)
