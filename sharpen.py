"""
Three ways to turn a coarse temperature map (1 km cells) into a fine one (10 m), using Sentinel-2
reflectance at 10 m as the guide:

    resample  enlarge the coarse cells smoothly (bilinear): no new information, the "do nothing" line
    tsharp    TsHARP: one straight line, LST = a + b x NDVI, fitted on the coarse cells, used at 10 m
    dms       data mining sharpener: an ensemble of regression trees on all ten bands plus NDVI, NDWI
              and NDMI, fitted on the most uniform coarse cells and used at 10 m

tsharp and dms end with residual correction: inside every coarse cell the fine values are shifted
so that their average is the measured coarse value again. The map gains detail but never disagrees
with the measurement it started from.

Two kinds of residual correction:
    constant  one shift per coarse cell: exact, but neighbouring cells get different shifts, which
              draws straight 1 km edges into the fine map
    smooth    the shifts are treated as values at the cell centres and blended bilinearly, then the
              small error this leaves in each cell's average is smoothed and added again, a few rounds,
              until the map is both free of edges and (almost) exactly on the measurement
"""
import numpy as np
from scipy import ndimage
from sklearn.ensemble import RandomForestRegressor

FACTOR = 100                                   # 10 m pixels along one side of a 1 km cell


def block_mean(a, f):
    """Mean over f x f blocks, ignoring NaN; and the share of each block that is not NaN."""
    h, w = a.shape[-2] // f, a.shape[-1] // f
    b = a[..., :h * f, :w * f].reshape(a.shape[:-2] + (h, f, w, f))
    ok = np.isfinite(b)
    n = ok.sum(axis=(-3, -1))
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(ok, b, 0).sum(axis=(-3, -1)) / n
    return np.where(n > 0, mean, np.nan).astype("float32"), (n / (f * f)).astype("float32")


def expand(c, f):
    """Each coarse cell repeated as an f x f block."""
    return np.repeat(np.repeat(c, f, axis=-2), f, axis=-1)


def features(refl):
    """(13, H, W) float32: ten reflectances (0-1) plus NDVI, NDWI (green vs NIR) and NDMI (NIR vs SWIR)."""
    r = refl.astype("float32") / 10000
    blue, green, red, nir, swir1 = r[0], r[1], r[2], r[3], r[8]
    with np.errstate(invalid="ignore", divide="ignore"):
        idx = [(nir - red) / (nir + red), (green - nir) / (green + nir), (nir - swir1) / (nir + swir1)]
    return np.concatenate([r, np.stack(idx)]).astype("float32")


def resample(coarse, f=FACTOR):
    """Bilinear enlargement through the cell centres; gaps filled from the nearest cell first."""
    filled = coarse.copy()
    if np.isnan(filled).any():
        idx = ndimage.distance_transform_edt(np.isnan(filled), return_distances=False, return_indices=True)
        filled = filled[tuple(idx)]
    h, w = coarse.shape
    rows = (np.arange(h * f) + 0.5) / f - 0.5
    cols = (np.arange(w * f) + 0.5) / f - 0.5
    out = np.empty((h * f, w * f), dtype="float32")
    for i0 in range(0, h * f, 1000):                          # in bands, to keep memory low
        rr, cc = np.meshgrid(rows[i0:i0 + 1000], cols, indexing="ij")
        out[i0:i0 + 1000] = ndimage.map_coordinates(filled, [rr, cc], order=1, mode="nearest")
    out[np.isnan(expand(coarse, f))] = np.nan
    return out


def residual_correct(fine, coarse, f=FACTOR, smooth=False, rounds=8):
    """Shift the fine values in each coarse cell so that their mean equals the coarse value:
    one constant per cell, or (smooth) bilinear shifts refined over a few rounds."""
    if not smooth:
        mean, _ = block_mean(fine, f)
        return fine + expand(coarse - mean, f)
    out = fine.copy()
    for _ in range(rounds):
        mean, _ = block_mean(out, f)
        out = out + resample(coarse - mean, f)
    return out


def edge_ratio(fine, f=FACTOR):
    """Mean temperature jump between neighbouring pixels across coarse-cell borders, divided by the
    mean jump between neighbours inside cells. About 1: no visible grid. Above 1: the grid shows."""
    across, inside = [], []
    for d in (np.abs(np.diff(fine, axis=1)), np.abs(np.diff(fine, axis=0)).T):
        border = (np.arange(d.shape[1]) % f) == f - 1          # pair (f-1, f) straddles a border
        across.append(d[:, border])
        inside.append(d[:, ~border])
    a = np.nanmean(np.concatenate([x.ravel() for x in across]))
    b = np.nanmean(np.concatenate([x.ravel() for x in inside]))
    return float(a / b)


def conservation_error(fine, coarse, f=FACTOR):
    """Mean absolute gap between each coarse cell's measured value and the mean of its fine pixels."""
    mean, _ = block_mean(fine, f)
    return float(np.nanmean(np.abs(mean - coarse)))


def _training_cells(coarse, feats, clear, f, keep):
    """Coarse cells to learn from: measured, at least 90% clear in Sentinel-2, and among the `keep`
    share that is most uniform inside (lowest spread of NDVI), so each cell's average is meaningful."""
    fine = np.where(clear, feats, np.nan)
    cf, share = block_mean(fine, f)
    ndvi = fine[10]
    h, w = coarse.shape
    spread = np.nanstd(ndvi[:h * f, :w * f].reshape(h, f, w, f), axis=(1, 3))
    ok = np.isfinite(coarse) & (share[0] >= 0.9) & np.isfinite(cf).all(0)
    cut = np.quantile(spread[ok], keep) if ok.any() else 0
    return ok & (spread <= cut), cf


def tsharp(coarse, feats, clear, f=FACTOR, keep=0.8):
    use, cf = _training_cells(coarse, feats, clear, f, keep)
    b, a = np.polyfit(cf[10][use], coarse[use], 1)
    fine = np.where(clear, a + b * feats[10], np.nan).astype("float32")
    return residual_correct(fine, coarse, f), {"a": float(a), "b": float(b), "cells": int(use.sum())}


def dms(coarse, feats, clear, f=FACTOR, keep=0.8, seed=0, smooth=False):
    raw, info = dms_raw(coarse, feats, clear, f, keep, seed)
    return residual_correct(raw, coarse, f, smooth=smooth), info


def dms_raw(coarse, feats, clear, f=FACTOR, keep=0.8, seed=0):
    """The trees' 10 m prediction before any residual correction."""
    use, cf = _training_cells(coarse, feats, clear, f, keep)
    model = RandomForestRegressor(n_estimators=60, min_samples_leaf=3, max_features=0.5, n_jobs=-1, random_state=seed)
    model.fit(cf[:, use].T, coarse[use])
    flat = feats.reshape(feats.shape[0], -1)
    good = clear.ravel() & np.isfinite(flat).all(0)
    pred = np.full(flat.shape[1], np.nan, dtype="float32")
    where = np.flatnonzero(good)
    for i in range(0, len(where), 1_000_000):                 # predict in chunks
        part = where[i:i + 1_000_000]
        pred[part] = model.predict(flat[:, part].T)
    fine = pred.reshape(clear.shape)
    imp = dict(zip(["B02", "B03", "B04", "B08", "B05", "B06", "B07", "B8A", "B11", "B12", "NDVI", "NDWI", "NDMI"],
                   np.round(model.feature_importances_, 3)))
    return fine, {"cells": int(use.sum()), "importance": imp}
