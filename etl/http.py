"""
HTTP helpers built on requests with sensible retries.
"""
from __future__ import annotations
import os
from pathlib import Path
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


def session_with_retries(total: int = 5, backoff: float = 0.3) -> requests.Session:
    s = requests.Session()
    retry = Retry(
        total=total,
        read=total,
        connect=total,
        backoff_factor=backoff,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=["GET", "HEAD"],
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    s.headers.update({"User-Agent": "sustai-geo-app/1.0"})
    return s


def head_probe(url: str, timeout: int = 20) -> requests.Response:
    """HEAD request to check content length/type before full download."""
    s = session_with_retries()
    r = s.head(url, timeout=timeout, allow_redirects=True)
    # Some servers don't support HEAD; fall back to GET with stream
    if r.status_code >= 400:
        r = s.get(url, timeout=timeout, stream=True)
    r.raise_for_status()
    return r


def stream_download(url: str, dst: str | Path, chunk_size: int = 1 << 20) -> Path:
    """Stream a large file to disk with retries and backoff."""
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)

    s = session_with_retries()
    with s.get(url, stream=True, timeout=60) as resp:
        resp.raise_for_status()
        with open(dst, "wb") as f:
            for chunk in resp.iter_content(chunk_size):
                if chunk:
                    f.write(chunk)

    # basic sanity
    if dst.stat().st_size == 0:
        raise IOError(f"Downloaded 0 bytes from {url}")
    return dst
