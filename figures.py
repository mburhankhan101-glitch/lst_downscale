"""
A figure that explains constant vs smooth residual correction, for the README.

    python figures.py 2026-05-06        # writes runs/correction_explained.png (a few minutes)

Top left: the idea with made-up numbers (two 1 km cells, measured 36 and 33 deg C, a model that says
34 in both). Top right: the edge ratio of both corrections, from runs/summary.csv. Bottom: real data
from one date, at the 1 km border in the study area where the two corrections differ most: the 10 m
map with each correction, and the average temperature walking across that border.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import lst_data
from sharpen import FACTOR, block_mean, dms_raw, features, residual_correct

OUT = "runs"


def idea(ax):
    x = np.arange(200)                                  # 200 pixels: 100 in cell A, 100 in cell B
    measured = np.where(x < 100, 36.0, 33.0)
    model = np.full(200, 34.0)
    constant = model + np.where(x < 100, 2.0, -1.0)
    centres = np.array([50, 150])
    smooth = model + np.interp(x, centres, [2.0, -1.0])
    ax.plot(x, model, color="grey", lw=1.5, label="model before correction (34)")
    ax.plot(x, constant, color="#c0392b", lw=2.5, label="constant correction: a cliff")
    ax.plot(x, smooth, color="#1f6aa3", lw=2.5, label="smooth correction: a ramp")
    for lo, hi, v in ((0, 100, 36), (100, 200, 33)):
        ax.hlines(v, lo, hi, colors="black", linestyles="dashed", lw=1)
    ax.text(50, 36.25, "cell A measured 36", ha="center", fontsize=9)
    ax.text(150, 33.25, "cell B measured 33", ha="center", fontsize=9)
    ax.axvline(100, color="black", lw=0.8, alpha=0.5)
    ax.text(101, 34.9, "1 km border", fontsize=8, rotation=90, va="bottom")
    ax.annotate("3 deg C jump between\ntwo pixels 10 m apart", xy=(100, 34.5), xytext=(122, 35.6), fontsize=8,
                color="#c0392b", arrowprops=dict(arrowstyle="->", color="#c0392b"))
    ax.set_ylim(32.3, 37)
    ax.set_xlabel("pixels across two 1 km cells (10 m each)")
    ax.set_ylabel("temperature (deg C)")
    ax.set_title("The idea: one shift per cell makes a cliff, a blended shift makes a ramp", fontsize=10)
    ax.legend(fontsize=8, loc="lower left")


def edge_bars(ax):
    s = pd.read_csv(os.path.join(OUT, "summary.csv"))
    s = s[s.test.str.startswith("B") & s.method.str.startswith("DMS")]
    vals = [float(s[s.method == m].edge_ratio.iloc[0]) for m in ("DMS trees", "DMS trees, smooth")]
    bars = ax.bar(["constant", "smooth"], vals, color=["#c0392b", "#1f6aa3"], width=0.55)
    ax.axhline(1.0, color="black", lw=0.8, ls="dashed")
    ax.text(0.5, 1.03, "1.0 = no seams", fontsize=8, va="bottom", ha="center")
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.03, f"{v:.2f}", ha="center", fontsize=10)
    ax.set_ylim(0, max(vals) + 0.3)
    ax.set_ylabel("jumps across 1 km borders / jumps inside cells")
    ax.set_title("Edge ratio, real Sentinel-3, mean of 5 dates", fontsize=10)


def real(date, axes):
    ls = lst_data.landsat(date)
    s3 = lst_data.sentinel3(date, ls["time"])
    s2 = lst_data.sentinel2(date)
    feats = features(s2["refl"])
    clear = s2["clear"] & np.isfinite(feats).all(0)
    raw, _ = dms_raw(s3["lst"], feats, clear)
    const = residual_correct(raw, s3["lst"])
    smooth = residual_correct(raw, s3["lst"], smooth=True)
    # the left-right pair of fully clear cells whose constant shifts differ most
    shift = s3["lst"] - block_mean(raw, FACTOR)[0]
    share = block_mean(np.where(clear, 1.0, np.nan), FACTOR)[1]
    jump = np.abs(np.diff(shift, axis=1))
    jump[(share[:, :-1] < 0.95) | (share[:, 1:] < 0.95)] = 0
    r, c = np.unravel_index(np.nanargmax(jump), jump.shape)
    rows = slice(r * FACTOR, (r + 1) * FACTOR)
    cols = slice(c * FACTOR, (c + 2) * FACTOR)           # two cells side by side, border in the middle
    lo, hi = np.nanpercentile(np.r_[const[rows, cols].ravel(), smooth[rows, cols].ravel()], [2, 98])
    for ax, img, name in ((axes[0], const, "constant"), (axes[1], smooth, "smooth")):
        ax.imshow(img[rows, cols], cmap="inferno", vmin=lo, vmax=hi, interpolation="nearest",
                  extent=(0, 2, 1, 0))
        ax.axvline(1, color="white", lw=0.8, ls="dashed")
        ax.set_title(f"Real 10 m map, {name} correction\n(2 x 1 km, dashed line = 1 km border)", fontsize=10)
        ax.set_xticks([]); ax.set_yticks([])
    x = (np.arange(2 * FACTOR) + 0.5) * 10 - FACTOR * 10
    for img, colour, name in ((const, "#c0392b", "constant"), (smooth, "#1f6aa3", "smooth")):
        axes[2].plot(x, np.nanmean(img[rows, cols], axis=0), color=colour, lw=2, label=name)
    axes[2].axvline(0, color="black", lw=0.8, alpha=0.5)
    axes[2].set_xlabel("metres from the 1 km border")
    axes[2].set_ylabel("mean temperature of each column (deg C)")
    axes[2].set_title(f"Walking across that border ({date})", fontsize=10)
    axes[2].legend(fontsize=8)


def main():
    date = sys.argv[1] if len(sys.argv) > 1 else "2026-05-06"
    fig = plt.figure(figsize=(15, 9.5), layout="constrained")
    grid = fig.add_gridspec(2, 3)
    idea(fig.add_subplot(grid[0, :2]))
    edge_bars(fig.add_subplot(grid[0, 2]))
    real(date, [fig.add_subplot(grid[1, k]) for k in range(3)])
    fig.suptitle("Constant vs smooth residual correction", fontsize=13)
    fig.savefig(os.path.join(OUT, "correction_explained.png"), dpi=80)
    print("Saved runs/correction_explained.png")


if __name__ == "__main__":
    main()
