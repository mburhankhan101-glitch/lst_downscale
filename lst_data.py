"""
The inputs for one date over the study area, on three aligned grids (UTM zone 43N):

    10 m   Sentinel-2 L2A reflectance, 10 bands (the 20 m bands resampled), clouds and shadows from SCL
    100 m  Landsat 8/9 Collection 2 surface temperature (thermal band measured at 100 m, delivered on a
           30 m grid), averaged onto 100 m cells; clouds, shadows and their margins from the QA band
    1 km   Sentinel-3 SLSTR land surface temperature: swath pixels (each with its own latitude and
           longitude) placed on 1 km cells; cloudy pixels and flagged pixels left out. Clouds come from
           the basic cloud tests (cloud_in); the Bayesian cloud flag marked 85% of a clear May day as
           cloudy over hot, bright ground, so it is not used, and the clear-sky masks of Landsat and
           Sentinel-2 are applied as well

The study area is a 40 x 40 km square around the crop_project fields south of Lahore (74.30 E,
31.19 N), its corners on whole kilometres so that 100 cells of 10 m make exactly one 1 km cell.
Everything comes from Microsoft Planetary Computer and is cached in cache/.

    python lst_data.py 2026-05-06       # download and cache one date, print what was found
"""
import os
import sys
from datetime import datetime, timedelta

import numpy as np
import planetary_computer
import pystac_client
import rasterio
import requests
import xarray as xr
from pyproj import Transformer
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject, transform_bounds
from rasterio.windows import Window
from scipy.interpolate import griddata
from scipy.spatial import cKDTree

CRS = "EPSG:32643"
CENTRE = (74.302, 31.193)                     # the crop_project fields
SIZE = 40_000                                 # metres
CACHE = "cache"
S2_BANDS = ["B02", "B03", "B04", "B08", "B05", "B06", "B07", "B8A", "B11", "B12"]
CLEAR_SCL = (4, 5, 6, 7)                      # vegetation, bare soil, water, unclassified
STAC = "https://planetarycomputer.microsoft.com/api/stac/v1"


def grids():
    """(x0, y_top) of the square and the transforms of the 10 m, 100 m and 1 km grids."""
    x, y = Transformer.from_crs("EPSG:4326", CRS, always_xy=True).transform(*CENTRE)
    x0 = np.floor((x - SIZE / 2) / 1000) * 1000
    top = np.floor((y + SIZE / 2) / 1000) * 1000
    return {res: (from_origin(x0, top, res, res), (SIZE // res, SIZE // res)) for res in (10, 100, 1000)}


def bbox_lonlat():
    tr, shape = grids()[1000]
    return transform_bounds(CRS, "EPSG:4326", *rasterio.transform.array_bounds(*shape, tr))


def catalog():
    return pystac_client.Client.open(STAC, modifier=planetary_computer.sign_inplace)


def _read(href, tr, shape, resampling, src_nodata=0):
    """One raster band onto our grid (only the window that covers it is downloaded)."""
    dst = np.zeros(shape, dtype="float32")
    with rasterio.open(href) as src:
        b = transform_bounds(CRS, src.crs, *rasterio.transform.array_bounds(*shape, tr))
        win = rasterio.windows.from_bounds(b[0] - 300, b[1] - 300, b[2] + 300, b[3] + 300, src.transform)
        try:
            win = win.intersection(Window(0, 0, src.width, src.height)).round_offsets().round_lengths()
        except rasterio.errors.WindowError:
            return dst
        reproject(src.read(1, window=win).astype("float32"), dst, src_transform=src.window_transform(win),
                  src_crs=src.crs, dst_transform=tr, dst_crs=CRS, resampling=resampling,
                  src_nodata=src_nodata, dst_nodata=0)
    return dst


def _cached(name, make):
    path = os.path.join(CACHE, name)
    if os.path.exists(path):
        return dict(np.load(path, allow_pickle=True))
    out = make()
    os.makedirs(CACHE, exist_ok=True)
    np.savez_compressed(path, **out)
    return out


def sentinel2(date):
    """Reflectance (10, H, W) uint16 x 10000 and a clear-sky mask on the 10 m grid, from the same
    day (or one day either side), clearest scenes first."""
    def make():
        tr, shape = grids()[10]
        d = datetime.fromisoformat(date)
        items = list(catalog().search(collections=["sentinel-2-l2a"], bbox=bbox_lonlat(),
                                      datetime=f"{(d - timedelta(days=1)).date()}/{(d + timedelta(days=1)).date()}").items())
        items.sort(key=lambda i: (abs((i.datetime.date() - d.date()).days), i.properties["eo:cloud_cover"]))
        refl = np.zeros((len(S2_BANDS),) + shape, dtype="uint16")
        scl = np.zeros(shape, dtype="uint8")
        used = []
        for item in items:
            empty = scl == 0
            if not empty.any():
                break
            offset = 1000 if item.properties.get("s2:processing_baseline", "00.00") >= "04.00" else 0
            s = _read(item.assets["SCL"].href, tr, shape, Resampling.nearest)
            take = empty & (s > 0)
            if not take.any():
                continue
            for k, band in enumerate(S2_BANDS):
                v = _read(item.assets[band].href, tr, shape, Resampling.bilinear)
                refl[k][take] = np.clip(v[take] - offset, 1, 60000).astype("uint16")
            scl[take] = s[take]
            used.append(item.id)
        return {"refl": refl, "clear": np.isin(scl, CLEAR_SCL), "items": np.array(used)}
    return _cached(f"s2_{date}.npz", make)


def landsat(date):
    """Surface temperature (deg C) and the share of clear 30 m pixels in each 100 m cell; the scene's time."""
    def make():
        tr, shape = grids()[100]
        items = list(catalog().search(collections=["landsat-c2-l2"], bbox=bbox_lonlat(), datetime=date,
                                      query={"platform": {"in": ["landsat-8", "landsat-9"]}}).items())
        items.sort(key=lambda i: i.properties["eo:cloud_cover"])
        lst = np.full(shape, np.nan, dtype="float32")
        share_bad = np.ones(shape, dtype="float32")
        used = []
        for item in items:                    # one pass is cut into scenes along the track: stitch them
            dn = _read(item.assets["lwir11"].href, tr, shape, Resampling.average)
            take = np.isnan(lst) & (dn > 0)
            if not take.any():
                continue
            lst[take] = dn[take] * 0.00341802 + 149.0 - 273.15
            # QA bits: 0 fill, 1 dilated cloud, 3 cloud, 4 cloud shadow, 5 snow -> not clear
            with rasterio.open(item.assets["qa_pixel"].href) as src:
                b = transform_bounds(CRS, src.crs, *rasterio.transform.array_bounds(*shape, tr))
                win = rasterio.windows.from_bounds(b[0] - 300, b[1] - 300, b[2] + 300, b[3] + 300, src.transform)
                win = win.intersection(Window(0, 0, src.width, src.height)).round_offsets().round_lengths()
                q = src.read(1, window=win)
                bad = ((q & 1) | ((q >> 1) & 1) | ((q >> 3) & 1) | ((q >> 4) & 1) | ((q >> 5) & 1)).astype("float32")
                part = np.ones(shape, dtype="float32")
                reproject(bad, part, src_transform=src.window_transform(win), src_crs=src.crs,
                          dst_transform=tr, dst_crs=CRS, resampling=Resampling.average)
            share_bad[take] = part[take]
            used.append(item.id)
        return {"lst": lst, "clear": np.isfinite(lst) & (share_bad < 0.01),
                "time": np.array(items[0].datetime.isoformat()), "item": np.array(" ".join(used))}
    return _cached(f"landsat_{date}.npz", make)


def _download(href, path):
    if not os.path.exists(path):
        r = requests.get(href, timeout=600)
        r.raise_for_status()
        with open(path, "wb") as f:
            f.write(r.content)
    return path


def sentinel3(date, near_time):
    """LST (deg C) on the 1 km grid from the daytime Sentinel-3 pass closest in time to near_time.

    Each swath pixel has its own latitude and longitude. A 1 km cell takes the value interpolated
    from the clear pixels around it, unless the swath pixel nearest to its centre is cloudy or flagged
    (so no cell borrows a value from across a cloud's edge)."""
    def make():
        tr, shape = grids()[1000]
        t = datetime.fromisoformat(str(near_time))
        items = [i for i in catalog().search(collections=["sentinel-3-slstr-lst-l2-netcdf"], bbox=bbox_lonlat(),
                                             datetime=date).items() if 3 <= i.datetime.hour <= 8]
        item = min(items, key=lambda i: abs((i.datetime - t).total_seconds()))
        os.makedirs(CACHE, exist_ok=True)
        files = {k: _download(item.assets[k].href, os.path.join(CACHE, f"{item.id[:31]}_{k}.nc"))
                 for k in ("lst-in", "slstr-geodetic-in", "slstr-flags-in")}
        lst = xr.open_dataset(files["lst-in"], engine="h5netcdf")
        geo = xr.open_dataset(files["slstr-geodetic-in"], engine="h5netcdf")
        # raw flags: these files declare 0 as their fill value, but 0 here means "no flag raised"
        flags = xr.open_dataset(files["slstr-flags-in"], engine="h5netcdf", mask_and_scale=False)
        exc = xr.open_dataset(files["lst-in"], engine="h5netcdf", mask_and_scale=False)["exception"].values
        lat, lon = geo["latitude_in"].values, geo["longitude_in"].values
        x, y = Transformer.from_crs("EPSG:4326", CRS, always_xy=True).transform(lon, lat)
        x0, top = tr.c, tr.f
        near = (x > x0 - 5000) & (x < x0 + SIZE + 5000) & (y < top + 5000) & (y > top - SIZE - 5000)
        k = lst["LST"].values
        good = (np.isfinite(k) & (flags["cloud_in"].values == 0) & (exc == 0))[near]
        px, py, val = x[near], y[near], (k[near] - 273.15)
        cx = x0 + (np.arange(shape[1]) + 0.5) * 1000
        cy = top - (np.arange(shape[0]) + 0.5) * 1000
        gx, gy = np.meshgrid(cx, cy)
        out = griddata((px[good], py[good]), val[good], (gx, gy), method="linear")
        _, nearest = cKDTree(np.c_[px, py]).query(np.c_[gx.ravel(), gy.ravel()])
        out[~good[nearest].reshape(shape)] = np.nan
        return {"lst": out.astype("float32"), "time": np.array(item.datetime.isoformat()), "item": np.array(item.id),
                "cloudy_share": np.array(1 - good.mean())}
    return _cached(f"s3_{date}.npz", make)


if __name__ == "__main__":
    for date in sys.argv[1:]:
        ls = landsat(date)
        s3 = sentinel3(date, ls["time"])
        s2 = sentinel2(date)
        print(f"{date}: Landsat {ls['item']} ({ls['clear'].mean():.0%} of 100 m cells clear), "
              f"Sentinel-3 {s3['item']} at {str(s3['time'])[11:16]} UTC ({np.isfinite(s3['lst']).mean():.0%} of 1 km cells), "
              f"Sentinel-2 {list(s2['items'])} ({s2['clear'].mean():.0%} clear)", flush=True)
