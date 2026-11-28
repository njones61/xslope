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
    # Pillow sizes its color table by the limit asked for, so ask small first:
    # a plot has a few thousand colors at most, and a full-size table for an
    # 8-frame sheet is a quarter gigabyte. The list is the same either way.
    counts = None
    for limit in (1 << 12, 1 << 16, sheet.size[0] * sheet.size[1]):
        counts = sheet.getcolors(min(limit, sheet.size[0] * sheet.size[1]))
        if counts is not None:
            break
    counts = counts or []
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


def _to_palette(im, pal, cache=None):
    """``im`` (RGB) as a P image on ``pal``, every pixel its nearest palette
    color. Done here rather than by Pillow's palette conversion, whose lookup
    is approximate (it turned white into 254, 254, 254).

    ``cache`` (a dict, kept across the frames of one file) remembers the
    nearest palette color of every RGB value already seen in a lookup table
    over all 2**24 colors, so a frame costs one table lookup per pixel and the
    nearest-color search only for colors no earlier frame had."""
    import numpy as np
    from PIL import Image
    if cache is None:
        cache = {}
    if "lut" not in cache:
        cache["lut"] = np.zeros(1 << 24, np.uint8)
        cache["known"] = np.zeros(1 << 24, bool)
    lut, known = cache["lut"], cache["known"]
    a = np.asarray(im, dtype=np.uint8).reshape(-1, 3)
    key = (a[:, 0].astype(np.int32) << 16) | (a[:, 1].astype(np.int32) << 8) | a[:, 2]
    new = np.unique(key[~known[key]])
    if new.size:
        cols = np.stack([(new >> 16) & 255, (new >> 8) & 255, new & 255], 1).astype(np.int32)
        p = pal.astype(np.int32)
        for s in range(0, len(new), 4096):
            d = ((cols[s:s + 4096, None, :] - p[None]) ** 2).sum(2)
            lut[new[s:s + 4096]] = np.argmin(d, axis=1)
        known[new] = True
    out = Image.fromarray(lut[key].reshape(im.size[1], im.size[0]), mode="P")
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

    The palette needs the sample, so frames are held until the end, but each is
    held as its RGB pixels compressed (a plot is mostly flat color) and is
    decompressed and put on the palette one at a time as the GIF is written, so
    no full-size copy of every frame is ever held (a 14 x 8 in, 200-frame trace
    peaked at 2.4 GB held the stock way). An exception while frames are being
    grabbed (a progress callback's Cancel) writes nothing.

    A GIF stores each frame's delay in hundredths of a second, so only frame
    rates that divide 100 (10, 20, 25, 50) play at exactly their rate; 15 fps
    plays at 60 ms a frame, 10% fast.
    """
    import contextlib
    import zlib
    from matplotlib import animation
    from PIL import Image

    class _Writer(animation.PillowWriter):
        def setup(self, fig, outfile, dpi=None):
            super().setup(fig, outfile, dpi=dpi)
            self._packed = []
            self._aborted = False

        def grab_frame(self, **savefig_kwargs):
            super().grab_frame(**savefig_kwargs)
            im = self._frames.pop().convert("RGB")
            self._packed.append((im.size, zlib.compress(im.tobytes(), 1)))

        @contextlib.contextmanager
        def saving(self, *args, **kwargs):
            with super().saving(*args, **kwargs) as w:
                try:
                    yield w
                except BaseException:
                    self._aborted = True
                    raise

        def finish(self):
            packed, self._packed = self._packed, []
            if self._aborted or not packed:
                return

            def unpack(item):
                size, data = item
                return Image.frombytes("RGB", size, zlib.decompress(data))

            pick = [unpack(it) for it in packed[:: max(1, len(packed) // 8)][:8]]
            w, h = pick[0].size
            sheet = Image.new("RGB", (w, h * len(pick)))
            for i, im in enumerate(pick):
                sheet.paste(im, (0, h * i))
            pal = _shared_palette(sheet, colors)
            del pick, sheet
            cache = {}
            first = _to_palette(unpack(packed[0]), pal, cache)
            rest = (_to_palette(unpack(it), pal, cache) for it in packed[1:])
            first.save(self.outfile, save_all=True, append_images=rest,
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
