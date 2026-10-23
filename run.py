"""
Sharpen 1 km land surface temperature to 10 m with three methods and check them against Landsat.

    Test A  Landsat (100 m) averaged to 1 km, sharpened back to 10 m, averaged to 100 m and compared
            with the Landsat it came from. The truth is known, so this measures the method alone.
    Test B  Real Sentinel-3 (1 km) sharpened to 10 m, averaged to 100 m and compared with Landsat of
            the same day. Adds the real-world gaps: two satellites, minutes apart, different angles and
            algorithms.

Methods: resample, TsHARP, and the tree sharpener (DMS) with constant or with smooth residual
correction (the same tree prediction, corrected two ways).

Two scores at 100 m, on cells clear in every source:
    RMSE / bias         how far each 100 m value is from Landsat (deg C)
    detail r, detail RMSE   the same, but only for the pattern INSIDE each 1 km cell (each value minus
                        its 1 km cell's mean). This is what sharpening is supposed to add; a method that
                        only copies the 1 km value scores r = 0 here.
Nothing can check the 10 m detail itself: Landsat's thermal band is 100 m. At 10 m two more numbers:
    edge ratio          jumps between neighbouring pixels across 1 km borders / jumps inside cells
                        (about 1: no visible grid)
    conservation error  mean gap between each 1 km measurement and the mean of its 10 m pixels

    python run.py 2026-03-27 2026-05-06 2026-05-14 2026-06-07 2026-06-23     # writes runs/
    python run.py --figures 2026-05-06                                         # redraw a figure only
    python run.py --tag heatwave2022 2022-05-27        # separate tables: runs/results_heatwave2022.csv ...
    python run.py --area lahore --tag lahore2022 2022-05-27      # Lahore district (see lst_data.AREAS)
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio

import lst_data
from sharpen import (FACTOR, block_mean, conservation_error, dms_raw, edge_ratio, expand, features, resample,
                     residual_correct, tsharp)

OUT = "runs"


def name(stem):
    """Output file name: the default area keeps the original names, others carry the area's name."""
    return stem if lst_data.AREA == "south" else f"{lst_data.AREA}_{stem}"


def scores(pred100, truth100, ok):
    """RMSE, bias and the within-1 km detail scores on the cells in ok."""
    p, t = np.where(ok, pred100, np.nan), np.where(ok, truth100, np.nan)
    d = p - t
    pm, _ = block_mean(p, 10)
    tm, _ = block_mean(t, 10)
    dp, dt = (p - expand(pm, 10))[ok], (t - expand(tm, 10))[ok]
    return {"rmse": float(np.sqrt(np.nanmean(d ** 2))), "bias": float(np.nanmean(d)),
            "detail_r": float(np.corrcoef(dp, dt)[0, 1]) if dp.std() > 0 else 0.0,
            "detail_rmse": float(np.sqrt(np.mean((dp - dt) ** 2))), "cells_100m": int(ok.sum())}


def run_date(date, s2=None):
    ls = lst_data.landsat(date)
    s3 = lst_data.sentinel3(date, ls["time"])
    s2 = lst_data.sentinel2(date) if s2 is None else s2       # another guide image can be passed in
    feats = features(s2["refl"])
    clear10 = s2["clear"] & np.isfinite(feats).all(0)
    truth = np.where(ls["clear"], ls["lst"], np.nan)
    _, s2_share = block_mean(np.where(clear10, 1.0, np.nan), 10)          # clear share of each 100 m cell
    coarse_a, share_a = block_mean(truth, 10)
    coarse_a[share_a < 0.9] = np.nan                                       # a fair 1 km average needs most of it
    tests = {"A: Landsat averaged to 1 km": coarse_a, "B: real Sentinel-3": s3["lst"]}
    rows, keep = [], {}
    for test, coarse in tests.items():
        raw, info_trees = dms_raw(coarse, feats, clear10)
        for name in ("resample", "TsHARP", "DMS trees", "DMS trees, smooth"):
            if name == "resample":
                fine, info = np.where(clear10, resample(coarse), np.nan), {}
            elif name == "TsHARP":
                fine, info = tsharp(coarse, feats, clear10)
            else:
                fine = residual_correct(raw, coarse, smooth=name.endswith("smooth"))
                info = info_trees
            p100, share = block_mean(fine, 10)
            ok = np.isfinite(truth) & np.isfinite(p100) & (share >= 0.5) & (s2_share >= 0.5) \
                & np.isfinite(expand(coarse_a, 10)) & np.isfinite(expand(s3["lst"], 10))
            rows.append({"date": date, "test": test, "method": name, **scores(p100, truth, ok),
                         "edge_ratio": edge_ratio(fine), "conservation_error": conservation_error(fine, coarse)})
            print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in rows[-1].items()}, flush=True)
            if name.startswith("DMS"):
                keep[(test, name)] = (fine, info)
    # Sentinel-3 against Landsat at 1 km: how different are the two sensors before any sharpening?
    both = np.isfinite(coarse_a) & np.isfinite(s3["lst"])
    diff = (s3["lst"] - coarse_a)[both]
    sensor = {"date": date, "s3_minus_landsat_1km_bias": float(diff.mean()),
              "s3_minus_landsat_1km_rmse": float(np.sqrt((diff ** 2).mean())),
              "s3_landsat_1km_r": float(np.corrcoef(s3["lst"][both], coarse_a[both])[0, 1]),
              "cells_1km": int(both.sum()), "s3_time": str(s3["time"])[11:16], "landsat_time": str(ls["time"])[11:16],
              "importance": json.dumps(keep[("B: real Sentinel-3", "DMS trees")][1].get("importance", {}))}
    return rows, sensor, (s2, ls, s3, truth, keep)


def figure(date, s2, ls, s3, truth, keep, write_tif=True):
    tr, _ = lst_data.grids()[10]
    fine_b = keep[("B: real Sentinel-3", "DMS trees, smooth")][0]
    fine_c = keep.get(("B: real Sentinel-3", "DMS trees"), (None,))[0]
    rgb = np.stack([s2["refl"][2], s2["refl"][1], s2["refl"][0]], -1).astype("float32") / 3000
    lo, hi = np.nanpercentile(truth, [2, 98])
    _, _, width, height = lst_data.extent_m()
    w_km, h_km = width // 1000, height // 1000
    area = lst_data.AREAS[lst_data.AREA]
    fig, ax = plt.subplots(2, 4, figsize=(19, 10), layout="constrained")
    rgb100 = np.stack([block_mean(rgb[..., k], 10)[0] for k in range(3)], -1)
    panels = [(rgb100.clip(0, 1), None, "Sentinel-2 true colour"),
              (s3["lst"], "inferno", "Sentinel-3 LST, 1 km (input)"),
              (truth, "inferno", "Landsat LST, 100 m (the check)"),
              (block_mean(fine_b, 10)[0], "inferno", "Sharpened Sentinel-3, averaged to 100 m")]
    for a, (img, cmap, title) in zip(ax.ravel()[:4], panels):
        a.imshow(img, cmap=cmap, vmin=None if cmap is None else lo, vmax=None if cmap is None else hi,
                 extent=(0, w_km, 0, h_km), interpolation="nearest")
        a.set_title(title, fontsize=10)
        a.set_xticks([]); a.set_yticks([])
    # zoom: 3 x 3 km around the area's point of interest at full 10 m resolution
    from pyproj import Transformer
    x, y = Transformer.from_crs("EPSG:4326", lst_data.CRS, always_xy=True).transform(*area["zoom"])
    col, row = ~tr * (x, y)
    r0, c0 = int(row) - 150, int(col) - 150
    z = (slice(r0, r0 + 300), slice(c0, c0 + 300))
    ax[1, 0].imshow(rgb[z].clip(0, 1), interpolation="nearest")
    ax[1, 0].set_title(f"{area['zoom_name']}, 3 x 3 km, true colour", fontsize=10)
    if fine_c is not None:
        ax[1, 1].imshow(fine_c[z], cmap="inferno", vmin=lo, vmax=hi, interpolation="nearest")
        ax[1, 1].set_title("10 m, constant correction\n(1 km block edges)", fontsize=10)
    im = ax[1, 2].imshow(fine_b[z], cmap="inferno", vmin=lo, vmax=hi, interpolation="nearest")
    ax[1, 2].set_title("10 m, smooth correction\n(detail inferred, not measured)", fontsize=10)
    for a in ax[1]:
        a.set_xticks([]); a.set_yticks([])
    ax[1, 3].axis("off")
    fig.colorbar(im, ax=ax, shrink=0.6, label="land surface temperature (deg C)")
    place = "South of Lahore" if lst_data.AREA == "south" else lst_data.AREA.capitalize()
    fig.suptitle(f"{place}, {date}, {w_km} x {h_km} km", fontsize=13)
    fig.savefig(os.path.join(OUT, name(f"lst_{date}.png")), dpi=75)
    plt.close(fig)
    if not write_tif:
        return
    with rasterio.open(os.path.join(OUT, name(f"lst_10m_{date}.tif")), "w", driver="GTiff", height=fine_b.shape[0],
                       width=fine_b.shape[1], count=1, dtype="float32", crs=lst_data.CRS, transform=tr,
                       nodata=np.nan, compress="deflate", predictor=3) as dst:
        dst.write(fine_b.astype("float32"), 1)


def main():
    os.makedirs(OUT, exist_ok=True)
    argv = sys.argv[1:]
    if argv[0] == "--area":
        lst_data.set_area(argv[1])
        argv = argv[2:]
    if argv[0] == "--figures":                     # redraw from the saved 10 m maps, no recomputing
        for date in argv[1:]:
            ls = lst_data.landsat(date)
            with rasterio.open(os.path.join(OUT, name(f"lst_10m_{date}.tif"))) as src:
                fine = src.read(1)
            figure(date, lst_data.sentinel2(date), ls, lst_data.sentinel3(date, ls["time"]),
                   np.where(ls["clear"], ls["lst"], np.nan), {("B: real Sentinel-3", "DMS trees, smooth"): (fine, {})}, write_tif=False)
        return
    args = argv
    tag = ""
    if args[0] == "--tag":
        tag, args = "_" + args[1], args[2:]
    rows, sensors = [], []
    for date in args:
        r, s, extra = run_date(date)
        rows += r
        sensors.append(s)
        figure(date, *extra)
        pd.DataFrame(rows).round(3).to_csv(os.path.join(OUT, f"results{tag}.csv"), index=False)
        pd.DataFrame(sensors).round(3).to_csv(os.path.join(OUT, f"sensors{tag}.csv"), index=False)
    res = pd.DataFrame(rows)
    summary = res.groupby(["test", "method"])[["rmse", "bias", "detail_r", "detail_rmse", "edge_ratio",
                                                "conservation_error"]].mean().round(2)
    summary.to_csv(os.path.join(OUT, f"summary{tag}.csv"))
    print(summary.to_string())
    print(pd.DataFrame(sensors).drop(columns="importance").round(2).to_string(index=False))


if __name__ == "__main__":
    main()
