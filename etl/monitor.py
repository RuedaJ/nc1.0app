"""
Lightweight resource monitoring helpers using psutil.
"""
from __future__ import annotations
import time
import contextlib
import psutil

_PROC = psutil.Process()


def mem_used_gb() -> float:
    return _PROC.memory_info().rss / (1024 ** 3)


@contextlib.contextmanager
def stage(name: str, warn_at_gb: float = 6.0):
    t0 = time.time()
    m0 = mem_used_gb()
    try:
        yield
    finally:
        dt = time.time() - t0
        m1 = mem_used_gb()
        delta = m1 - m0
        msg = f"[STAGE] {name}: {dt:.1f}s | ΔRSS≈{delta:.2f} GB | RSS≈{m1:.2f} GB"
        if m1 > warn_at_gb:
            msg = "[WARN] " + msg
        print(msg)
