"""
Simulate and plot a ChIP-seq-like signal track over a small genomic interval.

Produces a single-track "coverage/pileup" plot: noisy background reads plus
one (or more) enriched peak(s) with a realistic asymmetric shape, an optional
peak-call bar underneath (like a MACS2-style called-peak annotation), and a
summit marker.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# ---- palette (single categorical hue, light mode) ----
BLUE       = "#2a78d6"   # series / fill
BLUE_FILL  = "#2a78d6"
INK        = "#0b0b0b"   # primary text
INK_MUTED  = "#898781"   # axis / tick labels
BASELINE   = "#c3c2b7"   # axis line
GRID       = "#e1e0d9"   # hairline gridlines
SURFACE    = "#fcfcfb"


def simulate_chipseq_track(
    start=1_000_000,
    length=3_000,          # bp span of the interval
    n_points=1_500,        # resolution of the simulated track
    peak_center_frac=0.45, # peak summit as a fraction of the interval
    peak_width=180,        # bp, controls peak breadth
    peak_height=48,        # summit signal above background
    background_mean=3.0,   # background read depth (Poisson mean)
    skew=0.35,             # asymmetry of the peak (0 = symmetric)
    seed=7,
):
    """Return (positions, signal, peak_call_bounds, summit_pos)."""
    rng = np.random.default_rng(seed)

    x = np.linspace(start, start + length, n_points)
    center = start + peak_center_frac * length

    # Skewed-Gaussian peak shape: two half-widths (left/right) for asymmetry,
    # mimicking real ChIP-seq peaks which are rarely perfectly symmetric.
    sigma_left = peak_width * (1 - skew)
    sigma_right = peak_width * (1 + skew)
    sigma = np.where(x < center, sigma_left, sigma_right)
    peak_shape = np.exp(-0.5 * ((x - center) / sigma) ** 2)

    signal_mean = background_mean + peak_height * peak_shape

    # Poisson noise on top of the mean trace gives a realistic "jagged" pileup
    # look instead of a perfectly smooth curve.
    raw_counts = rng.poisson(lam=np.clip(signal_mean, 0.05, None)).astype(float)

    # Light smoothing (moving average) approximates the fragment-extension /
    # binning smoothing seen in genome browser tracks.
    window = 15
    kernel = np.ones(window) / window
    smoothed = np.convolve(raw_counts, kernel, mode="same")

    # Called-peak region: where the smoothed signal exceeds an enrichment
    # threshold over background, widened slightly like a real peak call.
    threshold = background_mean + 0.35 * peak_height
    above = smoothed > threshold
    if above.any():
        idx = np.where(above)[0]
        call_start, call_end = x[idx[0]], x[idx[-1]]
    else:
        call_start, call_end = center - peak_width, center + peak_width

    summit_pos = x[np.argmax(smoothed)]

    return x, smoothed, (call_start, call_end), summit_pos


def plot_chipseq_track(x, y, call_bounds, summit_pos, title="Simulated ChIP-seq signal"):
    fig, (ax_track, ax_call) = plt.subplots(
        nrows=2, figsize=(8, 3.6), sharex=True,
        gridspec_kw={"height_ratios": [5, 0.6], "hspace": 0.06},
    )
    fig.patch.set_facecolor(SURFACE)
    for ax in (ax_track, ax_call):
        ax.set_facecolor(SURFACE)

    # --- signal track: thin line + soft fill to baseline ---
    ax_track.fill_between(x, y, 0, color=BLUE_FILL, alpha=0.18, linewidth=0)
    ax_track.plot(x, y, color=BLUE, linewidth=1.6, solid_joinstyle="round")

    # summit marker
    summit_y = y[np.argmin(np.abs(x - summit_pos))]
    ax_track.plot([summit_pos], [summit_y], marker="o", markersize=5,
                  color=BLUE, markeredgecolor=SURFACE, markeredgewidth=1)
    ax_track.annotate(
        "summit", xy=(summit_pos, summit_y), xytext=(0, 8),
        textcoords="offset points", ha="center",
        fontsize=9, color=INK_MUTED,
    )

    ax_track.set_ylabel("Signal\n(reads per bin)", fontsize=9, color=INK_MUTED)
    ax_track.set_title(title, fontsize=11, color=INK, loc="left", pad=10)
    ax_track.spines[["top", "right"]].set_visible(False)
    ax_track.spines["left"].set_color(BASELINE)
    ax_track.spines["bottom"].set_visible(False)
    ax_track.tick_params(axis="y", colors=INK_MUTED, labelsize=8, length=0)
    ax_track.tick_params(axis="x", length=0)
    ax_track.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
    ax_track.set_ylim(bottom=0)
    ax_track.margins(x=0)

    # --- called-peak bar (rounded rectangle) beneath the track ---
    call_start, call_end = call_bounds
    ax_call.add_patch(FancyBboxPatch(
        (call_start, 0.15), call_end - call_start, 0.7,
        boxstyle="round,pad=0,rounding_size=25",
        linewidth=0, facecolor=BLUE,
    ))
    ax_call.set_ylim(0, 1)
    ax_call.set_yticks([])
    ax_call.set_ylabel("Peak\ncall", fontsize=9, color=INK_MUTED, rotation=0,
                        ha="right", va="center", labelpad=28)
    ax_call.spines[:].set_visible(False)
    ax_call.margins(x=0)

    # x-axis: genomic coordinates, muted ticks
    ax_call.set_xlabel("Genomic position (bp)", fontsize=9, color=INK_MUTED)
    ax_call.tick_params(axis="x", colors=INK_MUTED, labelsize=8, length=3, color=BASELINE)
    ax_call.xaxis.set_major_formatter(lambda v, pos: f"{v:,.0f}")

    fig.tight_layout()
    return fig


if __name__ == "__main__":
    x, y, call_bounds, summit_pos = simulate_chipseq_track()
    fig = plot_chipseq_track(x, y, call_bounds, summit_pos)
    fig.savefig("chipseq_peak_simulation.png", dpi=200, bbox_inches="tight")
    print("Saved chipseq_peak_simulation.png")
