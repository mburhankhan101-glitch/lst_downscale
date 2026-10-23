"""
Does the Sentinel-2 guide image have to be from the same day as Sentinel-3?

Sentinel-3 measures temperature every day, Sentinel-2 passes every 5 days and is often cloudy, so a
same-day pair is rare. A composite (the median of the clear observations over several days) has no
clouds, but it shows the land as it was over a period, not on the day: in May fields are harvested,
irrigated and ploughed within days. This experiment sharpens the SAME Sentinel-3 image (27 May 2022)
with different guides and scores each with Tests A and B, exactly as in run.py:

    same day                      the Sentinel-2 image of 27 May (the best case)
    one image, 5 days earlier     22 May
    composite, +-10 days          17 May to 6 June, 4 images
    composite, +-30 days          27 April to 26 June, 12 images
    composite, April              the month before, 6 images

Composites use only Sentinel-2 passes of the same relative orbit as 27 May (orbit 5), so the viewing
angle is the same and only time changes, and they never include 27 May itself. Two areas: the farmland
south of Lahore ("south", 40 x 40 km) and central Lahore around Data Darbar ("lahore_city", 20 x 20 km).
Also reported: how different each guide is from the same-day image (mean NDVI difference).

    python composites.py south
    python composites.py lahore_city
    python composites.py --figure            # runs/composites.png from both tables
"""
import os
import sys
from datetime import date as Date, timedelta

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import lst_data
import run
from sharpen import features

DATE, ORBIT = "2022-05-27", 5
OUT = "runs"


def guides():
    d = Date.fromisoformat(DATE)
    around = lambda a, b: lst_data.sentinel2_days(str(d + timedelta(days=a)), str(d + timedelta(days=b)),
                                                  ORBIT, skip=(DATE,))
    return [("same day", None),
            ("one image, 5 days earlier", ["2022-05-22"]),
            ("composite, +-10 days", around(-10, 10)),
            ("composite, +-30 days", around(-30, 30)),
            ("composite, April", lst_data.sentinel2_days("2022-04-01", "2022-04-30", ORBIT, skip=(DATE,)))]


def ndvi(s2):
    f = features(s2["refl"])[10]
    return np.where(s2["clear"], f, np.nan)


def experiment(area):
    lst_data.set_area(area)
    same = lst_data.sentinel2(DATE)
    ndvi_same = ndvi(same)
    rows = []
    for label, days in guides():
        if days is None:
            s2 = same
        elif len(days) == 1:
            s2 = lst_data.sentinel2(days[0], exact=True, orbit=ORBIT)
        else:
            s2 = lst_data.composite(days, ORBIT, f"{DATE}_{label.split(', ')[1].replace('+-', 'pm').replace(' ', '')}")
        both = np.isfinite(ndvi_same) & s2["clear"]
        diff = np.abs(ndvi(s2) - ndvi_same)[both]
        change = {"guide": label, "days": " ".join(days) if days else DATE, "images": len(days) if days else 1,
                  "guide_clear_share": float(s2["clear"].mean()),
                  "ndvi_diff_mean": float(diff.mean()) if diff.size else np.nan,
                  "ndvi_changed_share": float((diff > 0.2).mean()) if diff.size else np.nan}
        print(change, flush=True)
        r, _, _ = run.run_date(DATE, s2=s2)
        rows += [{"area": area, **change, **x} for x in r]
        pd.DataFrame(rows).round(3).to_csv(os.path.join(OUT, f"composites_{area}.csv"), index=False)
    return pd.DataFrame(rows)


def figure():
    df = pd.concat([pd.read_csv(os.path.join(OUT, f"composites_{a}.csv")) for a in ("south", "lahore_city")])
    df = df[df.method == "DMS trees, smooth"]
    order = [g for g, _ in [("same day", 0), ("one image, 5 days earlier", 0), ("composite, +-10 days", 0),
                            ("composite, +-30 days", 0), ("composite, April", 0)]]
    names = {"south": "farmland south of Lahore", "lahore_city": "central Lahore (Data Darbar)"}
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.5), layout="constrained")
    x = np.arange(len(order))
    for k, (col, title) in enumerate((("detail_r", "Pattern inside each 1 km cell explained (detail r)"),
                                      ("detail_rmse", "Error of that pattern (detail RMSE, deg C)"),
                                      ("ndvi_diff_mean", "How different the guide is from the same day\n(mean NDVI difference)"))):
        for j, area in enumerate(("south", "lahore_city")):
            for t, (test, hatch) in enumerate((("A: Landsat averaged to 1 km", ""), ("B: real Sentinel-3", "//"))):
                if col == "ndvi_diff_mean" and t == 1:
                    continue
                sub = df[(df.area == area) & (df.test == test)].set_index("guide").reindex(order)
                width = 0.2 if col != "ndvi_diff_mean" else 0.35
                off = (j * 2 + t - 1.5) * 0.2 if col != "ndvi_diff_mean" else (j - 0.5) * 0.35
                ax[k].bar(x + off, sub[col].values, width, color=["#2e7d4f", "#7b5ea7"][j], hatch=hatch,
                          edgecolor="white", label=f"{names[area]}" + ("" if col == "ndvi_diff_mean" else f", test {test[0]}"))
        ax[k].set_xticks(x, [o.replace("composite, ", "composite\n").replace("one image, ", "one image\n") for o in order],
                         fontsize=8)
        ax[k].set_title(title, fontsize=10)
        ax[k].grid(axis="y", alpha=0.3)
    ax[0].legend(fontsize=7, loc="lower left")
    ax[2].legend(fontsize=7)
    fig.suptitle(f"Sharpening Sentinel-3 of {DATE} with Sentinel-2 guides of different dates (trees, smooth correction)",
                 fontsize=12)
    fig.savefig(os.path.join(OUT, "composites.png"), dpi=90)
    print("Saved runs/composites.png")


if __name__ == "__main__":
    if sys.argv[1] == "--figure":
        figure()
    else:
        print(experiment(sys.argv[1]).groupby(["guide", "test", "method"])[["detail_r", "detail_rmse"]].mean().round(2))
