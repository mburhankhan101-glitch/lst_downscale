"""
Find the days when Sentinel-3 (LST), Sentinel-2 and Landsat 8/9 passed closest together over an area.

For every day in the window: the daytime Sentinel-3 SLSTR LST pass nearest in time to Sentinel-2,
Sentinel-2's time and cloud cover, and Landsat's time, cloud cover and how much of the area its
scenes cover that day. Days are ranked: all three on the same day first, then by the time gap between
Sentinel-3 and Sentinel-2. Cloud cover is the whole scene's figure from the catalogue, so a scene
can be clear over the area even with some cloud elsewhere (and the other way round).

    python find_dates.py lahore 2022-03-20 2022-05-30
    python find_dates.py islamabad 2022-03-20 2022-05-30
"""
import os
import sys
from collections import defaultdict

import pandas as pd
from pyproj import Transformer
from shapely.geometry import box, shape
from shapely.ops import transform

import lst_data

PLACES = {"lahore": (74.302, 31.193), "islamabad": (73.05, 33.68)}     # 40 x 40 km around each


def area(centre, size=40_000):
    to_utm = Transformer.from_crs("EPSG:4326", lst_data.CRS, always_xy=True)          # both places are in UTM 43N
    x, y = to_utm.transform(*centre)
    back = Transformer.from_crs(lst_data.CRS, "EPSG:4326", always_xy=True).transform
    return transform(back, box(x - size / 2, y - size / 2, x + size / 2, y + size / 2))


def find(place, start, end):
    aoi = area(PLACES[place])
    cat = lst_data.catalog()
    window = f"{start}/{end}"
    s3, s2, ls = defaultdict(list), defaultdict(list), defaultdict(list)
    for i in cat.search(collections=["sentinel-3-slstr-lst-l2-netcdf"], intersects=aoi, datetime=window).items():
        if 3 <= i.datetime.hour <= 8 and shape(i.geometry).contains(aoi):         # daytime, whole area in swath
            s3[i.datetime.date()].append(i)
    for i in cat.search(collections=["sentinel-2-l2a"], intersects=aoi, datetime=window).items():
        s2[i.datetime.date()].append(i)
    for i in cat.search(collections=["landsat-c2-l2"], intersects=aoi, datetime=window,
                        query={"platform": {"in": ["landsat-8", "landsat-9"]}}).items():
        ls[i.datetime.date()].append(i)
    rows = []
    for day in sorted(set(s3) & set(s2)):
        t2 = min(i.datetime for i in s2[day])
        best = min(s3[day], key=lambda i: abs((i.datetime - t2).total_seconds()))
        cover = shape(s2[day][0].geometry)
        for i in s2[day][1:]:
            cover = cover.union(shape(i.geometry))
        row = {"date": day, "s3_time": best.datetime.strftime("%H:%M"), "s3_satellite": best.id[:3],
               "s2_time": t2.strftime("%H:%M"), "s3_s2_gap_min": round(abs((best.datetime - t2).total_seconds()) / 60),
               "s2_cloud_pct": round(min(i.properties["eo:cloud_cover"] for i in s2[day]), 1),
               "s2_area_covered_pct": round(100 * cover.intersection(aoi).area / aoi.area)}
        if day in ls:
            lcover = shape(ls[day][0].geometry)
            for i in ls[day][1:]:
                lcover = lcover.union(shape(i.geometry))
            tl = ls[day][0].datetime
            row.update({"landsat": " ".join(sorted({i.properties["platform"].replace("landsat-", "L") for i in ls[day]})),
                        "landsat_time": tl.strftime("%H:%M"),
                        "s3_landsat_gap_min": round(abs((best.datetime - tl).total_seconds()) / 60),
                        "landsat_cloud_pct": round(max(i.properties["eo:cloud_cover"] for i in ls[day]), 1),
                        "landsat_area_covered_pct": round(100 * lcover.intersection(aoi).area / aoi.area)})
        rows.append(row)
    df = pd.DataFrame(rows)
    df["all_three"] = df.get("landsat_area_covered_pct", pd.Series(dtype=float)).fillna(0) >= 99
    return df.sort_values(["all_three", "s3_s2_gap_min", "s2_cloud_pct"], ascending=[False, True, True])


if __name__ == "__main__":
    place, start, end = sys.argv[1], sys.argv[2], sys.argv[3]
    df = find(place, start, end)
    os.makedirs("runs", exist_ok=True)
    out = os.path.join("runs", f"dates_{place}_{start}_{end}.csv")
    df.to_csv(out, index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(df.to_string(index=False))
    print(f"Saved {out}")
