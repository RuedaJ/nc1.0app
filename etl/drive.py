"""
Download utilities for Google Drive and generic HTTP URLs.

- parse_drive_url: extract file id from a Drive link.
- download_drive_any: best-effort download for Drive (uses gdown if present).
- download_to: generic entry that uses gdown for Drive or requests streaming otherwise.
"""

from __future__ import annotations
import re
from pathlib import Path
from typing import Optional

try:
    import gdown  # robust Google Drive downloader
except Exception:
    gdown = None

from etl.http import head_probe, stream_download


_DRIVE_FILE_PATTERNS = [
    re.compile(r"https?://drive\.google\.com/file/d/([^/]+)/view"),
    re.compile(r"https?://drive\.google\.com/uc\?id=([^&]+)"),
    re.compile(r"https?://drive\.google\.com/open\?id=([^&]+)"),
]


def parse_drive_url(url: str) -> Optional[str]:
    """Return Google Drive file id if the URL looks like a Drive link, else None."""
    for pat in _DRIVE_FILE_PATTERNS:
        m = pat.search(url)
        if m:
            return m.group(1)
    return None


def download_drive_any(url: str, dst: str | Path) -> Path:
    """
    Download a file from Google Drive using gdown if available.
    Raises if gdown is not available.
    """
    if gdown is None:
        raise RuntimeError("gdown not available; cannot download from Google Drive.")
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    # gdown supports fuzzy parsing of Drive urls
    out = gdown.download(url=url, output=str(dst), quiet=False, fuzzy=True)
    if out is None:
        raise IOError(f"gdown failed to download: {url}")
    return Path(out)


def download_to(url: str, dst: str | Path) -> Path:
    """
    Download any HTTP/Drive URL to dst.
    Uses gdown when link is a Google Drive URL and gdown is available,
    otherwise falls back to robust requests streaming.
    """
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)

    # Sanity-check: HEAD probe (non-fatal if it fails)
    try:
        resp = head_probe(url)
        size = int(resp.headers.get("Content-Length", "0") or 0)
        if size and size > 20 * (1 << 30):  # 20 GB guard
            raise ValueError(f"Refusing to download huge file ({size/1e9:.1f} GB).")
    except Exception as e:
        print(f"[drive] HEAD probe failed ({e}); continuing to download...")

    # Drive?
    if parse_drive_url(url) and gdown is not None:
        return download_drive_any(url, dst)

    # Generic HTTP
    return stream_download(url, dst)
