"""
Sharpen 1 km land surface temperature to 10 m with three methods and check them against Landsat.

    Test A  Landsat (100 m) averaged to 1 km, sharpened back to 10 m, averaged to 100 m and compared
            with the Landsat it came from. The truth is known, so this measures the method alone.
    Test B  Real Sentinel-3 (1 km) sharpened to 10 m, averaged to 100 m and compared with Landsat of
            the same day. Adds the real-world gaps: two satellites, minutes apart, different angles and
            algorithms.

Two scores at 100 m, on cells clear in every source:
    RMSE / bias         how far each 100 m value is from Landsat (deg C)
    detail r, detail RMSE   the same, but only for the pattern INSIDE each 1 km cell (each value minus
                        its 1 km cell's mean). This is what sharpening is supposed to add; a method that
                        only copies the 1 km value scores r = 0 here.
Nothing can check the 10 m detail itself: Landsat's thermal band is 100 m.

    python run.py 2026-03-27 2026-05-06 2026-05-14 2026-06-07 2026-06-23     # writes runs/
    python run.py --figures 2026-05-06                                         # redraw a figure only
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
from sharpen import FACTOR, block_mean, dms, expand, features, resample, tsharp

OUT = "runs"
FIELDS = (74.302, 31.193)


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


def run_date(date):
    ls = lst_data.landsat(date)
    s3 = lst_data.sentinel3(date, ls["time"])
    s2 = lst_data.sentinel2(date)
    feats = features(s2["refl"])
    clear10 = s2["clear"] & np.isfinite(feats).all(0)
    truth = np.where(ls["clear"], ls["lst"], np.nan)
    _, s2_share = block_mean(np.where(clear10, 1.0, np.nan), 10)          # clear share of each 100 m cell
    coarse_a, share_a = block_mean(truth, 10)
    coarse_a[share_a < 0.9] = np.nan                                       # a fair 1 km average needs most of it
    tests = {"A: Landsat averaged to 1 km": coarse_a, "B: real Sentinel-3": s3["lst"]}
    rows, keep = [], {}
    for test, coarse in tests.items():
        for name, method in (("resample", None), ("TsHARP", tsharp), ("DMS trees", dms)):
            if method is None:
                fine, info = np.where(clear10, resample(coarse), np.nan), {}
            else:
                fine, info = method(coarse, feats, clear10)
            p100, share = block_mean(fine, 10)
            ok = np.isfinite(truth) & np.isfinite(p100) & (share >= 0.5) & (s2_share >= 0.5) \
                & np.isfinite(expand(coarse_a, 10)) & np.isfinite(expand(s3["lst"], 10))
            rows.append({"date": date, "test": test, "method": name, **scores(p100, truth, ok)})
            print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in rows[-1].items()}, flush=True)
            if name == "DMS trees":
                keep[test] = (fine, info)
    # Sentinel-3 against Landsat at 1 km: how different are the two sensors before any sharpening?
    both = np.isfinite(coarse_a) & np.isfinite(s3["lst"])
    diff = (s3["lst"] - coarse_a)[both]
    sensor = {"date": date, "s3_minus_landsat_1km_bias": float(diff.mean()),
              "s3_minus_landsat_1km_rmse": float(np.sqrt((diff ** 2).mean())),
              "s3_landsat_1km_r": float(np.corrcoef(s3["lst"][both], coarse_a[both])[0, 1]),
              "cells_1km": int(both.sum()), "s3_time": str(s3["time"])[11:16], "landsat_time": str(ls["time"])[11:16],
              "importance": json.dumps(keep["B: real Sentinel-3"][1].get("importance", {}))}
    return rows, sensor, (s2, ls, s3, truth, keep)


def figure(date, s2, ls, s3, truth, keep, write_tif=True):
    tr, _ = lst_data.grids()[10]
    fine_b = keep["B: real Sentinel-3"][0]
    rgb = np.stack([s2["refl"][2], s2["refl"][1], s2["refl"][0]], -1).astype("float32") / 3000
    lo, hi = np.nanpercentile(truth, [2, 98])
    fig, ax = plt.subplots(2, 3, figsize=(15, 10.5), layout="constrained")
    rgb100 = np.stack([block_mean(rgb[..., k], 10)[0] for k in range(3)], -1)
    panels = [(rgb100.clip(0, 1), None, "Sentinel-2 true colour"),
              (s3["lst"], "inferno", "Sentinel-3 LST, 1 km (input)"),
              (truth, "inferno", "Landsat LST, 100 m (the check)"),
              (block_mean(fine_b, 10)[0], "inferno", "Sharpened Sentinel-3, averaged to 100 m")]
    for a, (img, cmap, title) in zip(ax.ravel()[:4], panels):
        a.imshow(img, cmap=cmap, vmin=None if cmap is None else lo, vmax=None if cmap is None else hi,
                 extent=(0, 40, 0, 40), interpolation="nearest")
        a.set_title(title, fontsize=10)
        a.set_xticks([]); a.set_yticks([])
    # zoom: 3 x 3 km around the crop_project fields at full 10 m resolution
    from pyproj import Transformer
    x, y = Transformer.from_crs("EPSG:4326", lst_data.CRS, always_xy=True).transform(*FIELDS)
    col, row = ~tr * (x, y)
    r0, c0 = int(row) - 150, int(col) - 150
    z = (slice(r0, r0 + 300), slice(c0, c0 + 300))
    ax[1, 1].imshow(rgb[z].clip(0, 1), interpolation="nearest")
    ax[1, 1].set_title("crop_project fields, 3 x 3 km, true colour", fontsize=10)
    im = ax[1, 2].imshow(fine_b[z], cmap="inferno", vmin=lo, vmax=hi, interpolation="nearest")
    ax[1, 2].set_title("Same 3 x 3 km, sharpened LST at 10 m\n(detail inferred, not measured)", fontsize=10)
    for a in (ax[1, 1], ax[1, 2]):
        a.set_xticks([]); a.set_yticks([])
    fig.colorbar(im, ax=ax, shrink=0.6, label="land surface temperature (deg C)")
    fig.suptitle(f"South of Lahore, {date}, 40 x 40 km", fontsize=13)
    fig.savefig(os.path.join(OUT, f"lst_{date}.png"), dpi=75)
    plt.close(fig)
    if not write_tif:
        return
    with rasterio.open(os.path.join(OUT, f"lst_10m_{date}.tif"), "w", driver="GTiff", height=fine_b.shape[0],
                       width=fine_b.shape[1], count=1, dtype="float32", crs=lst_data.CRS, transform=tr,
                       nodata=np.nan, compress="deflate", predictor=3) as dst:
        dst.write(fine_b.astype("float32"), 1)


def main():
    os.makedirs(OUT, exist_ok=True)
    if sys.argv[1] == "--figures":                 # redraw from the saved 10 m maps, no recomputing
        for date in sys.argv[2:]:
            ls = lst_data.landsat(date)
            with rasterio.open(os.path.join(OUT, f"lst_10m_{date}.tif")) as src:
                fine = src.read(1)
            figure(date, lst_data.sentinel2(date), ls, lst_data.sentinel3(date, ls["time"]),
                   np.where(ls["clear"], ls["lst"], np.nan), {"B: real Sentinel-3": (fine, {})}, write_tif=False)
        return
    rows, sensors = [], []
    for date in sys.argv[1:]:
        r, s, extra = run_date(date)
        rows += r
        sensors.append(s)
        figure(date, *extra)
        pd.DataFrame(rows).round(3).to_csv(os.path.join(OUT, "results.csv"), index=False)
        pd.DataFrame(sensors).round(3).to_csv(os.path.join(OUT, "sensors.csv"), index=False)
    res = pd.DataFrame(rows)
    summary = res.groupby(["test", "method"])[["rmse", "bias", "detail_r", "detail_rmse"]].mean().round(2)
    summary.to_csv(os.path.join(OUT, "summary.csv"))
    print(summary.to_string())
    print(pd.DataFrame(sensors).drop(columns="importance").round(2).to_string(index=False))


if __name__ == "__main__":
    main()
