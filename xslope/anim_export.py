"""One animation exporter for every XSLOPE animation.

GIF is always available (matplotlib's PillowWriter). MP4 needs ffmpeg; when
ffmpeg is not on the system, asking for one raises a plain error naming the
cause, and :func:`ffmpeg_available` lets a caller (Studio's save dialog) offer
MP4 only when it can be written.

Both the steady flow trace (:func:`xslope.plot_seep.animate_flow_trace`) and the
transient frame export (:func:`xslope.plot_seep.save_transient_animation`) write
through :func:`save_animation`.
"""

from pathlib import Path

FORMATS = (".gif", ".mp4")


def _shared_palette(sheet, colors):
    """A shared palette for ``sheet`` as an (n, 3) uint8 array: its most frequent
    colors kept exactly (the background, the fills, a legend swatch), the rest
    of the slots a median cut of everything, so flat areas keep their color."""
    import numpy as np
    from PIL import Image
    counts = sheet.getcolors(sheet.size[0] * sheet.size[1]) or []
    counts.sort(reverse=True)
    n_exact = min(max(colors - 16, colors // 2), len(counts))
    exact = [c for _, c in counts[:n_exact]]
    n_rest = colors - len(exact)
    rest = []
    if n_rest > 0:
        q = sheet.quantize(colors=n_rest, method=Image.Quantize.MEDIANCUT)
        flat = q.getpalette()[: 3 * n_rest]
        rest = [tuple(flat[i:i + 3]) for i in range(0, len(flat), 3)]
    return np.array(exact + rest, dtype=np.uint8)


def _to_palette(im, pal):
    """``im`` (RGB) as a P image on ``pal``, every pixel its nearest palette
    color. Done here rather than by Pillow's palette conversion, whose lookup
    is approximate (it turned white into 254, 254, 254)."""
    import numpy as np
    from PIL import Image
    a = np.asarray(im, dtype=np.uint8).reshape(-1, 3)
    key = (a[:, 0].astype(np.int32) << 16) | (a[:, 1].astype(np.int32) << 8) | a[:, 2]
    uniq, inv = np.unique(key, return_inverse=True)
    cols = np.stack([(uniq >> 16) & 255, (uniq >> 8) & 255, uniq & 255], 1).astype(np.int32)
    best = np.empty(len(uniq), np.uint8)
    p = pal.astype(np.int32)
    for s in range(0, len(uniq), 4096):
        d = ((cols[s:s + 4096, None, :] - p[None]) ** 2).sum(2)
        best[s:s + 4096] = np.argmin(d, axis=1)
    out = Image.fromarray(best[inv].reshape(im.size[1], im.size[0]), mode="P")
    flat = pal.ravel().tolist()
    out.putpalette(flat + [0] * (768 - len(flat)))
    return out


def _compact_gif_writer(fps, colors):
    """A PillowWriter that puts every frame on ONE shared palette, undithered.

    The stock writer lets Pillow pick an adaptive palette per frame and dither
    to it; dither noise differs from frame to frame and defeats the GIF's LZW
    compression. One palette built from a sample of the frames (their most
    frequent colors exactly, then a median cut), applied without dithering,
    keeps flat fills flat and true (a 200-frame flow trace went from 16.5 MB to
    10.5 MB); anti-aliased edges snap to the nearest of ``colors`` shades.

    A GIF stores each frame's delay in hundredths of a second, so only frame
    rates that divide 100 (10, 20, 25, 50) play at exactly their rate; 15 fps
    plays at 60 ms a frame, 10% fast.
    """
    from matplotlib import animation
    from PIL import Image

    class _Writer(animation.PillowWriter):
        def finish(self):
            frames = [im.convert("RGB") for im in self._frames]
            if not frames:
                return
            pick = frames[:: max(1, len(frames) // 8)][:8]
            w, h = pick[0].size
            sheet = Image.new("RGB", (w, h * len(pick)))
            for i, im in enumerate(pick):
                sheet.paste(im, (0, h * i))
            pal = _shared_palette(sheet, colors)
            out = [_to_palette(im, pal) for im in frames]
            out[0].save(self.outfile, save_all=True, append_images=out[1:],
                        duration=int(round(1000 / self.fps)), loop=0, optimize=False)

    return _Writer(fps=fps)


def ffmpeg_available():
    """True when matplotlib can find an ffmpeg executable to write MP4."""
    from matplotlib import animation
    try:
        return bool(animation.writers.is_available("ffmpeg"))
    except Exception:
        return False


def available_formats():
    """The file extensions this system can write: always .gif, .mp4 with ffmpeg."""
    return FORMATS if ffmpeg_available() else (".gif",)


def save_animation(anim, path, fps, dpi=None, progress=None, gif_colors=64):
    """Write a matplotlib animation to ``path``; the extension picks the format.

    Parameters
    ----------
    anim : matplotlib.animation.Animation
    path : str or Path ending in .gif or .mp4
    fps : frames per second of the written file
    dpi : resolution (None keeps the figure's own)
    progress : optional callable(frame_index, total_frames)
    gif_colors : the GIF's shared palette size (2-256)

    Returns the path written. Raises ValueError for another extension and
    RuntimeError for .mp4 when ffmpeg is not installed.
    """
    from matplotlib import animation
    path = Path(path)
    ext = path.suffix.lower()
    if ext == ".gif":
        writer = _compact_gif_writer(fps, int(min(max(gif_colors, 2), 256)))
    elif ext == ".mp4":
        if not ffmpeg_available():
            raise RuntimeError(
                "MP4 needs ffmpeg, which was not found on this system; "
                "save as .gif instead (or install ffmpeg).")
        writer = animation.FFMpegWriter(
            fps=fps, codec="h264", extra_args=["-pix_fmt", "yuv420p"])
    else:
        raise ValueError(f"animation format must be one of {FORMATS}, got {ext!r}")
    path.parent.mkdir(parents=True, exist_ok=True)
    kw = {"writer": writer}
    if dpi is not None:
        kw["dpi"] = dpi
    if progress is not None:
        kw["progress_callback"] = progress
    anim.save(str(path), **kw)
    return path
