"""Draw the background of the macOS disk image's Finder window.

The window shows the app on the left and an Applications shortcut on the right.
This background adds an arrow between them and two lines of instructions below:
drag the app onto Applications, then open it from there. The positions come
from the layout constants in ``build_app.py``, which also places the icons, so the
arrow always sits between them.

Writes ``dmg_background.png`` (1x) and ``dmg_background@2x.png`` (Retina) next to
this script; ``build_app.py`` combines the pair into one image for the disk image.
Run it on a Mac (it uses the system font) after changing the layout, and commit
both PNGs:

    python packaging/macos/make_dmg_background.py

A window with a background picture is always drawn in light mode, with black
icon labels, so a light background keeps the labels readable in dark mode too.
"""
import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from build_app import (APP_NAME, DMG_BACKGROUND, DMG_ICON_SIZE, DMG_ICONS,  # noqa: E402
                       DMG_LABEL_SIZE, DMG_WINDOW)

FONT = "/System/Library/Fonts/SFNS.ttf"
BACKGROUND = (255, 255, 255)
ARROW = (118, 118, 118)      # #767676: 4.5:1 on white
TEXT = (51, 51, 51)          # #333333: 12.6:1 on white
TEXT_SIZE = 18               # points
# Dragging only copies the app, and nothing launches, so the second line says
# what comes next.
INSTRUCTIONS = (f"To install, drag {APP_NAME} onto the Applications folder.",
                "Then open it from Applications or Launchpad.")


def draw(scale):
    w, h = DMG_WINDOW
    img = Image.new("RGB", (w * scale, h * scale), BACKGROUND)
    d = ImageDraw.Draw(img)
    (ax, ay), (bx, _by) = DMG_ICONS["app"], DMG_ICONS["Applications"]
    half = DMG_ICON_SIZE / 2

    # The arrow: from just past the app icon to just short of Applications.
    gap = 0.25 * DMG_ICON_SIZE
    x0, x1 = ax + half + gap, bx - half - gap
    shaft = 6
    head_len, head_half = 26, 18
    s = scale
    d.rectangle([x0 * s, (ay - shaft / 2) * s, (x1 - head_len) * s, (ay + shaft / 2) * s],
                fill=ARROW)
    d.polygon([((x1 - head_len) * s, (ay - head_half) * s), (x1 * s, ay * s),
               ((x1 - head_len) * s, (ay + head_half) * s)], fill=ARROW)

    # The instructions, starting one line below the icon labels. They are placed
    # from the labels, not from the window's bottom edge, which a path bar can cover.
    font = ImageFont.truetype(FONT, TEXT_SIZE * scale)
    label_bottom = ay + half + 2.2 * DMG_LABEL_SIZE
    for i, line in enumerate(INSTRUCTIONS):
        ty = label_bottom + (1.5 + 1.4 * i) * TEXT_SIZE
        d.text((w / 2 * s, ty * s), line, font=font, fill=TEXT, anchor="mm")
    return img


def main():
    one = DMG_BACKGROUND
    two = one.with_name(f"{one.stem}@2x{one.suffix}")
    for scale, path in ((1, one), (2, two)):
        img = draw(scale)
        img.save(path, dpi=(72 * scale, 72 * scale))
        print(f"wrote {path.relative_to(HERE.parent.parent)} ({img.width}x{img.height})")


if __name__ == "__main__":
    main()
