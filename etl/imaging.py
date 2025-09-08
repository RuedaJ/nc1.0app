"""
Small imaging utilities for report polishing using Pillow.
"""
from __future__ import annotations
from pathlib import Path
from typing import Optional
from PIL import Image, ImageDraw, ImageFont


def add_footer_badge(
    img_path: str | Path,
    out_path: str | Path,
    text: str = "sustai-geo-app",
    logo_path: Optional[str | Path] = None,
    opacity: int = 230,
) -> Path:
    img_path, out_path = Path(img_path), Path(out_path)
    im = Image.open(img_path).convert("RGBA")
    W, H = im.size
    badge_h = max(40, int(0.06 * H))

    canvas = Image.new("RGBA", (W, H + badge_h), (255, 255, 255, 0))
    canvas.paste(im, (0, 0))

    draw = ImageDraw.Draw(canvas)
    # footer bar
    draw.rectangle([(0, H), (W, H + badge_h)], fill=(255, 255, 255, opacity))

    # text
    try:
        font = ImageFont.truetype("DejaVuSans.ttf", size=max(14, badge_h // 2))
    except Exception:
        font = ImageFont.load_default()
    draw.text((12, H + (badge_h - font.size) // 2), text, fill=(0, 0, 0, 255), font=font)

    # logo (optional)
    if logo_path:
        logo = Image.open(logo_path).convert("RGBA")
        target = badge_h - 10
        logo.thumbnail((target, target))
        canvas.paste(logo, (W - logo.width - 12, H + (badge_h - logo.height) // 2), logo)

    out_path = out_path.with_suffix(".jpg")
    canvas.convert("RGB").save(out_path, quality=92)
    return out_path


def thumbnail(img_path: str | Path, out_path: str | Path, max_side: int = 512) -> Path:
    im = Image.open(img_path).convert("RGB")
    im.thumbnail((max_side, max_side))
    out = Path(out_path).with_suffix(".jpg")
    im.save(out, quality=88)
    return out
