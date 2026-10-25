"""
Do crops at the tail of a canal distributary run hotter, or green less, than crops at its head?

Distributaries (and minors) come from OpenStreetMap (data/osm_canals_south.json, (c) OpenStreetMap
contributors). Pieces with the same name are joined into one line; lines shorter than 15 km are
dropped. Lines keep their full length, so 0 and 1 are the real head and tail even where they lie
outside the area; only the parts inside the area have pixels. OSM draws waterways in the direction of flow; as a check, the head is taken to
be the end that lies closer to another canal (where a distributary leaves its parent canal), and it is
recorded whether that agrees with the drawing direction.

Every 20 m pixel within 1.5 km of a distributary (and nearer to it than to any other distributary) is
placed along it: 0 at the head, 1 at the tail. Pixels within 30 m of any canal and open water are left
out. For each line and date:

    green share         share of land pixels with NDVI >= 0.4 (crops), head third vs tail third
    green-crop anomaly  temperature of a green pixel minus the mean of the green pixels within 5 km of
                        it: how much hotter a crop is than crops around it, which removes the regional
                        gradient. Computed twice: with the downscaled 10 m map, and with Landsat's own
                        100 m measurement (independent of the downscaling model)
    tail - head         mean anomaly of the tail third minus the head third (deg C)

Placebo: 300 random straight lines of the same lengths, pointing in random directions (full circle),
scored the same way; and 300 more pointing the same way as the real distributaries (head to tail,
+-20 degrees), which also controls for a plain east-west difference across the area. Each distributary is compared with the spread of the placebo lines; and the mean
of all distributaries is compared with the means of as many random lines drawn 10,000 times (p: the
share of random sets whose mean is at least as far from zero), which tests whether the distributaries
together lean one way more than chance would.

    python canal.py               # writes runs/canal/ (about 25 minutes)
    python canal.py --figure      # redraw the figure from the saved tables
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
import shapely
from pyproj import Transformer
from scipy import ndimage
from shapely.geometry import LineString, MultiLineString
from shapely.ops import linemerge, transform

import lst_data
from sharpen import block_mean, expand

DATES = ["2026-03-27", "2026-05-06", "2026-05-14", "2026-06-07", "2026-06-23", "2022-05-27"]
OUT = os.path.join("runs", "canal")
STRIP, NEAR_CANAL, MIN_KM, STEP = 1500, 30, 15, 2          # metres, metres, km, pixels (2 x 10 m)
GREEN = 0.4


def distributaries():
    ways = json.load(open(os.path.join("data", "osm_canals_south.json")))
    to_utm = Transformer.from_crs("EPSG:4326", lst_data.CRS, always_xy=True).transform
    lines = [(w.get("tags", {}).get("name", ""), transform(to_utm, LineString([(p["lon"], p["lat"]) for p in w["geometry"]])))
             for w in ways]
    every = shapely.union_all([g for _, g in lines])
    out = []
    for name in sorted({n for n, _ in lines if any(k in n.lower() for k in ("distrib", "minor"))}):
        merged = linemerge(MultiLineString([g for n, g in lines if n == name]))
        line = max(getattr(merged, "geoms", [merged]), key=lambda g: g.length)
        if line.length < MIN_KM * 1000:
            continue
        others = shapely.union_all([g for n, g in lines if n != name])
        a, b = shapely.Point(line.coords[0]), shapely.Point(line.coords[-1])
        agrees = a.distance(others) <= b.distance(others)
        out.append({"name": name, "line": line if agrees else LineString(line.coords[::-1]),
                    "km": line.length / 1000, "osm_direction_agrees": agrees})
    return out, every


def neighbourhood_anomaly(lst, green, size_km=5):
    """Temperature minus the mean of green pixels within about size_km (computed on a 100 m grid)."""
    g100, _ = block_mean(np.where(green, lst, np.nan), 10)
    ok = np.isfinite(g100)
    n = int(size_km * 10) | 1
    s = ndimage.uniform_filter(np.where(ok, g100, 0), n, mode="constant")
    c = ndimage.uniform_filter(ok.astype("float32"), n, mode="constant")
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(c > 0.05, s / c, np.nan)
    return lst - expand(mean.astype("float32"), 10)


def pixel_positions(lines, shape, tr):
    """For sampled pixels near the lines: which line, fraction along it (0 head, 1 tail)."""
    rows, cols = np.mgrid[0:shape[0]:STEP, 0:shape[1]:STEP]
    xs, ys = tr * (cols + 0.5, rows + 0.5)
    pts = shapely.points(xs.ravel(), ys.ravel())
    dist = np.stack([shapely.distance(l, pts) for l in lines])
    which = dist.argmin(0)
    near = dist.min(0) <= STRIP
    frac = np.full(len(pts), np.nan)
    for k, l in enumerate(lines):
        sel = near & (which == k)
        frac[sel] = shapely.line_locate_point(l, pts[sel]) / l.length
    return rows.ravel(), cols.ravel(), which, near, frac, pts


def line_positions(line, tr, shape):
    """Sampled pixels within STRIP of one line (looking only around it): row, col, fraction along."""
    x0, y0, x1, y1 = line.buffer(STRIP).bounds
    d = 10 * STEP
    j0, j1 = max(int((x0 - tr.c) / d), 0), min(int((x1 - tr.c) / d) + 1, shape[1] // STEP)
    i0, i1 = max(int((tr.f - y1) / d), 0), min(int((tr.f - y0) / d) + 1, shape[0] // STEP)
    rows, cols = np.mgrid[i0:i1, j0:j1]
    xs, ys = tr.c + (cols * STEP + 0.5) * 10, tr.f - (rows * STEP + 0.5) * 10
    pts = shapely.points(xs.ravel(), ys.ravel())
    keep = shapely.distance(line, pts) <= STRIP
    frac = shapely.line_locate_point(line, pts[keep]) / line.length
    return rows.ravel()[keep] * STEP, cols.ravel()[keep] * STEP, frac


def score(frac, anom, green, land):
    """Tail third minus head third: green-crop anomaly (deg C) and green share."""
    head, tail = frac < 1 / 3, frac > 2 / 3
    a = lambda m: np.nanmean(anom[m & green]) if (m & green).sum() > 200 else np.nan
    g = lambda m: green[m & land].mean() if (m & land).sum() > 200 else np.nan
    return a(tail) - a(head), g(head), g(tail)


def random_lines(lengths, x0, top, width, height, n, rng, angles=None):
    """Random straight lines; with angles, each points like one of the real lines (head to tail, +-20 deg)."""
    out = []
    while len(out) < n:
        length = rng.choice(lengths)
        if angles is None:
            ang = rng.uniform(0, 2 * np.pi)             # the full circle: "tails" point every way
        else:
            ang = rng.choice(angles) + np.radians(rng.uniform(-20, 20))
        cx, cy = rng.uniform(x0, x0 + width), rng.uniform(top - height, top)
        dx, dy = np.cos(ang) * length / 2, np.sin(ang) * length / 2
        l = LineString([(cx - dx, cy - dy), (cx + dx, cy + dy)])
        if l.within(shapely.box(x0, top - height, x0 + width, top)):
            out.append(l)
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    tr, shape = lst_data.grids()[10]
    x0, top, width, height = lst_data.extent_m()
    lines, every = distributaries()
    print(f"{len(lines)} distributaries / minors of at least {MIN_KM} km:", flush=True)
    for d in lines:
        print(f"  {d['name']}: {d['km']:.1f} km, OSM direction {'agrees' if d['osm_direction_agrees'] else 'reversed'}", flush=True)
    rows, cols, which, near, frac, pts = pixel_positions([d["line"] for d in lines], shape, tr)
    off_canal = shapely.distance(every, pts) > NEAR_CANAL
    rng = np.random.default_rng(0)
    fakes = random_lines([d["line"].length for d in lines], x0, top, width, height, 300, rng)
    fake_pos = [line_positions(l, tr, shape) for l in fakes]
    heading = [np.arctan2(d["line"].coords[-1][1] - d["line"].coords[0][1], d["line"].coords[-1][0] - d["line"].coords[0][0])
               for d in lines]
    oriented = random_lines([d["line"].length for d in lines], x0, top, width, height, 300, rng, heading)
    oriented_pos = [line_positions(l, tr, shape) for l in oriented]

    summary, profiles, fakes_by = [], [], {}
    for date in DATES:
        s2 = lst_data.sentinel2(date)
        r = s2["refl"].astype("float32")
        with np.errstate(invalid="ignore", divide="ignore"):
            ndvi = (r[3] - r[2]) / (r[3] + r[2])
            ndwi = (r[1] - r[3]) / (r[1] + r[3])
        land = s2["clear"] & (ndwi < 0)
        green = land & (ndvi >= GREEN)
        with rasterio.open(os.path.join("runs", f"lst_10m_{date}.tif")) as src:
            lst10 = src.read(1)
        ls = lst_data.landsat(date)
        lst100 = expand(np.where(ls["clear"], ls["lst"], np.nan), 10)
        anoms = {"downscaled 10 m": neighbourhood_anomaly(lst10, green), "Landsat 100 m": neighbourhood_anomaly(lst100, green)}
        G, L = green[rows, cols] & off_canal, land[rows, cols] & off_canal
        for src_name, anom in anoms.items():
            A = anom[rows, cols]
            fake = np.array([score(fr, anom[r, c], green[r, c], land[r, c]) for r, c, fr in fake_pos])
            fakes_by[(date, src_name, "any direction")] = fake
            fakes_by[(date, src_name, "same direction as the canals")] = np.array(
                [score(fr, anom[r, c], green[r, c], land[r, c]) for r, c, fr in oriented_pos])
            p_lo, p_hi = np.nanpercentile(fake[:, 0], [5, 95])
            g_lo, g_hi = np.nanpercentile(fake[:, 2] - fake[:, 1], [5, 95])
            for k, d in enumerate(lines):
                sel = near & (which == k)
                diff, g_head, g_tail = score(frac[sel], A[sel], G[sel], L[sel])
                summary.append({"date": date, "temperature": src_name, "distributary": d["name"], "km": round(d["km"], 1),
                                "tail_minus_head_C": diff, "green_share_head": g_head, "green_share_tail": g_tail,
                                "placebo_5pct": p_lo, "placebo_95pct": p_hi,
                                "beyond_placebo": bool(np.isfinite(diff) and (diff > p_hi or diff < p_lo)),
                                "green_change": g_tail - g_head, "green_placebo_5pct": g_lo, "green_placebo_95pct": g_hi,
                                "green_beyond_placebo": bool(np.isfinite(g_tail - g_head) and (g_tail - g_head < g_lo or g_tail - g_head > g_hi))})
                bins = np.minimum((frac[sel] * 10).astype(int), 9)
                for b in range(10):
                    m = bins == b
                    profiles.append({"date": date, "temperature": src_name, "distributary": d["name"], "bin": b,
                                     "anomaly_C": np.nanmean(A[sel][m & G[sel]]) if (m & G[sel]).sum() > 50 else np.nan,
                                     "green_share": G[sel][m & L[sel]].mean() if (m & L[sel]).sum() > 50 else np.nan})
        print(f"{date}: done", flush=True)
    s = pd.DataFrame(summary)
    s.round(3).to_csv(os.path.join(OUT, "canal_summary.csv"), index=False)
    p = pd.DataFrame(profiles)
    p.round(3).to_csv(os.path.join(OUT, "canal_profiles.csv"), index=False)
    pd.DataFrame([{k: v for k, v in d.items() if k != "line"} for d in lines]).round(1).to_csv(
        os.path.join(OUT, "distributaries.csv"), index=False)
    figure(lines, every, p, s, tr, shape)
    by = s.groupby(["temperature", "date"]).agg(tail_minus_head_C=("tail_minus_head_C", "mean"),
                                                 lines_hotter_at_tail=("tail_minus_head_C", lambda x: int((x > 0).sum())),
                                                 lines=("tail_minus_head_C", lambda x: int(x.notna().sum())),
                                                 beyond_placebo=("beyond_placebo", "sum"),
                                                 placebo_5pct=("placebo_5pct", "first"), placebo_95pct=("placebo_95pct", "first"),
                                                 green_head=("green_share_head", "mean"), green_tail=("green_share_tail", "mean"),
                                                 lines_less_green_at_tail=("green_change", lambda x: int((x < 0).sum())),
                                                 green_beyond_placebo=("green_beyond_placebo", "sum"),
                                                 green_placebo_5pct=("green_placebo_5pct", "first"),
                                                 green_placebo_95pct=("green_placebo_95pct", "first"))
    tests = []
    for (date, src_name, placebo), fake in fakes_by.items():
        real = s[(s.date == date) & (s.temperature == src_name)]
        for col, values in (("tail_minus_head_C", fake[:, 0]), ("green_change", fake[:, 2] - fake[:, 1])):
            r = real[col].dropna().values
            v = values[np.isfinite(values)]
            if len(r) < 3:
                continue
            means = np.array([rng.choice(v, len(r), replace=False).mean() for _ in range(10_000)])
            tests.append({"date": date, "temperature": src_name, "measure": col, "placebo": placebo, "lines": len(r),
                          "real_mean": r.mean(), "random_mean_5pct": np.percentile(means, 5),
                          "random_mean_95pct": np.percentile(means, 95),
                          "p_two_sided": float(np.mean(np.abs(means - np.mean(v)) >= abs(r.mean() - np.mean(v))))})
    tests = pd.DataFrame(tests)
    tests.round(3).to_csv(os.path.join(OUT, "canal_tests.csv"), index=False)
    print(tests.round(3).to_string(index=False))
    by.round(2).to_csv(os.path.join(OUT, "canal_by_date.csv"))
    print(by.round(2).to_string())


def figure(lines, every, p, s, tr, shape):
    with rasterio.open(os.path.join("runs", "lst_10m_2022-05-27.tif")) as src:
        lst = src.read(1)
    x0, top, width, height = lst_data.extent_m()
    km = lambda g: ((np.array(g.xy[0]) - x0) / 1000, (np.array(g.xy[1]) - (top - height)) / 1000)
    fig, ax = plt.subplots(1, 3, figsize=(19, 6.3), layout="constrained")
    im = ax[0].imshow(block_mean(lst, 4)[0], cmap="inferno", extent=(0, width / 1000, 0, height / 1000),
                      vmin=np.nanpercentile(lst, 2), vmax=np.nanpercentile(lst, 98))
    for g in getattr(every, "geoms", [every]):
        ax[0].plot(*km(g), color="white", lw=0.4, alpha=0.5)
    for d in lines:
        xs, ys = km(d["line"])
        ax[0].plot(xs, ys, color="cyan", lw=1.6)
        ax[0].plot(xs[0], ys[0], "o", color="lime", ms=6)
        ax[0].plot(xs[-1], ys[-1], "s", color="red", ms=6)
    ax[0].set_xlim(-12, width / 1000 + 12)
    ax[0].set_ylim(-12, height / 1000 + 12)
    ax[0].set_xlabel("km"); ax[0].set_ylabel("km")
    ax[0].set_title("Distributaries (cyan) on the 27 May 2022 temperature map\ngreen dot = head, red square = tail", fontsize=10)
    fig.colorbar(im, ax=ax[0], shrink=0.7, label="LST (deg C)")
    for k, (col, title) in enumerate((("anomaly_C", "Green crops: hotter than crops within 5 km (deg C)"),
                                      ("green_share", "Share of land that is green crop"))):
        a = ax[k + 1]
        for date, colour in zip(DATES, plt.cm.viridis(np.linspace(0, 0.9, len(DATES)))):
            q = p[(p.date == date) & (p.temperature == "Landsat 100 m")].groupby("bin")[col].mean()
            a.plot((q.index + 0.5) / 10, q.values, marker="o", color=colour, label=date)
        a.axhline(0 if col == "anomaly_C" else np.nan, color="grey", lw=0.8)
        a.set_xlabel("position along the distributary (0 = head, 1 = tail)")
        a.set_title(title + "\n(mean of all distributaries; temperature from Landsat 100 m)", fontsize=10)
        a.legend(fontsize=8)
    fig.savefig(os.path.join(OUT, "canal.png"), dpi=85)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--figure":          # redraw from the saved tables only
        tr, shape = lst_data.grids()[10]
        lines, every = distributaries()
        figure(lines, every, pd.read_csv(os.path.join(OUT, "canal_profiles.csv")),
               pd.read_csv(os.path.join(OUT, "canal_summary.csv")), tr, shape)
    else:
        main()
