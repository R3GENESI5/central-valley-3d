"""Vapour-transport winds, column vapour and cloud base on storm days and on the afternoon after, winters 1959-60 to 1979-80
against 1999-2000 to 2021-22, for the 3D Central Valley page (Ali, 2026-11-18: "can we run the winds (that bring moisture
as model over it, i.e. like in earthnull school, showing TCWV, TPW, LCL, CAPE, before vs now" ... "ERA5 1959-80 vs 2000-22").
Page data, not a test. Reanalysis: a before/now difference mixes land change and weather; ERA5's land model has no irrigation,
but its 2 m temperature and dew point are analysed from station reports, so the floor's air carries what the stations saw.

Days (rule of src/corridor/corridor_ad.py, unchanged):
  day            the PRISM day, 12 UTC to 12 UTC (six-hourly steps 18, 00, 06, 12 UTC); only days with all four steps.
  inflow day     offshore inflow >= 150 kg/m/s: ERA5 column vapour flux along bearing 71 degrees, mean of the corridor's first
                 three 25 km bins (0-75 km from 37.80 N, 123.50 W, within 40 km of the line), mean of the four steps.
  storm days     all steps of inflow days, November to March.
  after storm    the 00 UTC step (16:00 PST, afternoon) of the first day without inflow after an inflow day, November to March.
  eras           water years 1960-1980 (winters 1959-60 to 1979-80) and 2000-2022 (winters 1999-2000 to 2021-22).
Fields:
  vapour transport  ERA5 viwve, viwvn (0.25 degree, era5_threshold), kg/m/s; the particles follow it.
  cloud base (LCL)  125 m per K of t2m - d2m (0.25 degree, era5_threshold), m above the ground.
  column vapour     ERA5 tcwv (1.5 degree, WeatherBench-2 era5_wb2_mimo2), kg/m2 = mm; TCWV and TPW are the same quantity.
  CAPE              not on disk; not shown.
Grid: the page's mesh at 5 km (x = x0 + 5000 i, y = y1 - 5000 j, EPSG:3310), bilinear from the ERA5 grid.
Writes figures/valley_3d/data/wind/wind.bin (float32; era, composite, field, row, column) and wind_meta.json.
Floor means (5 km cells at least half CVHM2 floor) and the offshore inflow, with water-year resample intervals, are printed.
Run  python src/landscape/valley3d_wind.py
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import xarray as xr
from PIL import Image
from rasterio.warp import transform
from scipy.interpolate import RegularGridInterpolator

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "src", "corridor"))
from track import BEARING, BIN_KM, HALF_KM, along_cross  # noqa: E402

Z = r"Z:\Data\california_central_valley"
DATA = os.path.join(REPO, "figures", "valley_3d", "data")
OUT = os.path.join(DATA, "wind")
ERAS = {"before": range(1960, 1981), "now": range(2000, 2023)}
COMPS = ["storm", "after"]
FIELDS = ["ivt_e", "ivt_n", "lcl", "tcwv"]
STEP = 5000.0
NBOOT, SEED = 2000, 20261118

meta = json.load(open(os.path.join(DATA, "meta.json"), encoding="utf-8"))
G = meta["grid"]
WX, WY = (G["mw"] - 1) // 5 + 1, (G["mh"] - 1) // 5 + 1
gx = G["x0"] + STEP * np.arange(WX)
gy = G["y1"] - STEP * np.arange(WY)
XX, YY = np.meshgrid(gx, gy)
lon, lat = transform("EPSG:3310", "EPSG:4326", XX.ravel().tolist(), YY.ravel().tolist())
LON, LAT = np.array(lon).reshape(WY, WX), np.array(lat).reshape(WY, WX)

# floor share of each 5 km cell, from floor.png (R + 256 G - 1 >= 0 is a CVHM2 floor cell, 300 m)
fl = np.asarray(Image.open(os.path.join(DATA, "floor.png")).convert("RGBA")).astype(np.int32)
isf = (fl[..., 0] + 256 * fl[..., 1] - 1) >= 0
k = int(STEP / G["tex_m"])
FLOOR = np.zeros((WY, WX), bool)
for j in range(WY):
    for i in range(WX):
        b = isf[max(0, int((j - 0.5) * k)):int((j + 0.5) * k), max(0, int((i - 0.5) * k)):int((i + 0.5) * k)]
        FLOOR[j, i] = b.size > 0 and b.mean() >= 0.5


def to_grid(field, la, lo):
    """field (lat, lon) on a regular grid; la, lo 1-D (lat may descend) -> page 5 km grid."""
    if la[0] > la[-1]:
        la, field = la[::-1], field[::-1]
    f = RegularGridInterpolator((la, lo), field, bounds_error=False, fill_value=np.nan)
    return f(np.c_[LAT.ravel(), LON.ravel()]).reshape(WY, WX)


def pday(t):
    return (pd.DatetimeIndex(t) + pd.Timedelta(hours=11, minutes=59)).normalize()


def main():
    if os.path.exists(os.path.join(OUT, "wind.bin")):
        raise SystemExit(f"{OUT}/wind.bin exists; not overwritten")
    os.makedirs(OUT, exist_ok=True)
    b = np.radians(BEARING)
    # per water year: sums and counts of each field per composite on the page grid, and offshore inflow and floor sums
    S = {e: {c: {f: np.zeros((WY, WX)) for f in FIELDS} for c in COMPS} for e in ERAS}
    N = {e: {c: 0 for c in COMPS} for e in ERAS}
    YR = []   # rows: era, wy, comp, n, floor lcl sum, floor tcwv sum, offshore inflow sum (per day for storm)
    for era, wys in ERAS.items():
        for wy in wys:
            parts = []
            for y in (wy - 1, wy):
                ds = xr.open_dataset(os.path.join(Z, "era5_threshold", f"era5_ivt_t2_td2_{y}.nc"))
                parts.append(ds.sel(time=ds.time.dt.month.isin([10, 11, 12] if y == wy - 1 else [1, 2, 3, 4])).load())
            T = xr.concat(parts, "time")
            mm = [xr.open_dataset(os.path.join(Z, "era5_wb2_mimo2", f"wb2_mimo2_{y}.nc"))[["tcwv"]] for y in (wy - 1, wy)]
            W = xr.concat([m.sel(time=m.time.dt.month.isin([10, 11, 12] if y == wy - 1 else [1, 2, 3, 4])).load() for m, y in zip(mm, (wy - 1, wy))], "time")
            la, lo = T.latitude.values, T.longitude.values
            LO2, LA2 = np.meshgrid(lo, la)
            al, cr = along_cross(LA2, LO2)
            kb = np.where((al >= 0) & (al < 3 * BIN_KM) & (np.abs(cr) <= HALF_KM), np.floor(al / BIN_KM).astype(int), -1)
            fa = T.viwve.values * np.sin(b) + T.viwvn.values * np.cos(b)
            f0 = np.mean([fa[:, kb == q].mean(1) for q in range(3)], axis=0)
            t = pd.DatetimeIndex(T.time.values)
            pd_ = pday(t)
            dfl = pd.DataFrame({"f": f0, "d": pd_})
            g = dfl.groupby("d").f.agg(["mean", "size"])
            g = g[g["size"] == 4]
            g = g[g.index.month.isin([11, 12, 1, 2, 3]) & (np.where(g.index.month >= 10, g.index.year + 1, g.index.year) == wy)]
            inflow = set(g.index[g["mean"] >= 150])
            after = set()
            for d in sorted(inflow):
                n1 = d + pd.Timedelta(days=1)
                if n1 in g.index and n1 not in inflow:
                    after.add(n1)
            sel = {"storm": np.array([d in inflow for d in pd_]),
                   "after": np.array([(d in after) and (h == 0) for d, h in zip(pd_, t.hour)])}
            wt = pd.DatetimeIndex(W.time.values)
            for c in COMPS:
                idx = np.where(sel[c])[0]
                if not len(idx):
                    YR.append((era, wy, c, 0, 0.0, 0.0, 0.0)); continue
                ts = t[idx]
                widx = [wt.get_loc(x) for x in ts]
                e_ = to_grid(T.viwve.values[idx].mean(0), la, lo)
                n_ = to_grid(T.viwvn.values[idx].mean(0), la, lo)
                l_ = to_grid((125.0 * (T.t2m.values[idx] - T.d2m.values[idx])).mean(0), la, lo)
                w_ = to_grid(W.tcwv.transpose("time", "latitude", "longitude").values[widx].mean(0), W.latitude.values, W.longitude.values - 360.0)
                for f, v in zip(FIELDS, (e_, n_, l_, w_)):
                    S[era][c][f] += v * len(idx)
                N[era][c] += len(idx)
                YR.append((era, wy, c, len(idx), np.nanmean(l_[FLOOR]) * len(idx), np.nanmean(w_[FLOOR]) * len(idx),
                           float(f0[idx].sum())))
            print(era, wy, "inflow days", len(inflow), "after-storm afternoons", len(after), flush=True)
    A = np.zeros((2, 2, 4, WY, WX), np.float32)
    for ei, e in enumerate(ERAS):
        for ci, c in enumerate(COMPS):
            for fi, f in enumerate(FIELDS):
                A[ei, ci, fi] = S[e][c][f] / max(1, N[e][c])
    A.tofile(os.path.join(OUT, "wind.bin"))
    Y = pd.DataFrame(YR, columns=["era", "wy", "comp", "n", "lcl", "tcwv", "f0"])
    Y.to_csv(os.path.join(OUT, "wind_years.csv"), index=False)
    rng = np.random.default_rng(SEED)
    res = {}
    for c in COMPS:
        for q in ("lcl", "tcwv", "f0"):
            y0, y1 = Y[(Y.era == "before") & (Y.comp == c)], Y[(Y.era == "now") & (Y.comp == c)]
            r = lambda d: d[q].sum() / d.n.sum()
            bs = [r(y1.iloc[rng.integers(0, len(y1), len(y1))]) - r(y0.iloc[rng.integers(0, len(y0), len(y0))]) for _ in range(NBOOT)]
            res[f"{c}_{q}"] = [round(r(y0), 1), round(r(y1), 1), round(r(y1) - r(y0), 1), round(float(np.quantile(bs, .025)), 1), round(float(np.quantile(bs, .975)), 1)]
            lab = {"lcl": "floor cloud base (m)", "tcwv": "floor column vapour (mm)", "f0": "offshore inflow (kg/m/s)"}[q]
            print(f"{c:5s} {lab:28s} before {res[f'{c}_{q}'][0]} | now {res[f'{c}_{q}'][1]} | change {res[f'{c}_{q}'][2]:+} [{res[f'{c}_{q}'][3]:+}, {res[f'{c}_{q}'][4]:+}]")
    steps = {e: {c: int(N[e][c]) for c in COMPS} for e in ERAS}
    print("six-hourly steps used:", steps, "; floor 5 km cells:", int(FLOOR.sum()))
    json.dump({"SOURCE": __doc__.strip(), "grid": {"wx": WX, "wy": WY, "step_m": STEP, "every_mesh": 5},
               "order": {"era": list(ERAS), "comp": COMPS, "field": FIELDS}, "steps": steps,
               "eras": {"before": "winters 1959-60 to 1979-80", "now": "winters 1999-2000 to 2021-22"},
               "floor": res, "floor_cells": int(FLOOR.sum())},
              open(os.path.join(OUT, "wind_meta.json"), "w", encoding="utf-8"), ensure_ascii=False)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
