"""Model wind and near-surface air over today's land and over the old land, same winter weather, for the 3D page
(Ali, 2026-11-18: "is there a way to show the wind now vs then, and the audience can see the temperature rise in low winds as they
pass over fallow, crops, and urban surfaces ... if it is modelled" ... "Yes, after drift read"). Page data, not a test.

Runs: the A28 WRF pair, winter 2013-14, d02 (12 km): ctrl (today's land) and pre_w (pre-1900 land), the pair read by the science
thread's three-ground profile read and by src/landscape/valley3d_lcl.py (same files, same floor-cell rule, same storm-hour rule).
Composites, 2013-11-01 to 2014-03-30 (150 days):
  afternoon   21 UTC (13 PST) of every day
  predawn     12 UTC (04 PST) of every day
  storm       hours with CTRL floor-mean rain >= 0.1 mm/h (rule of valley3d_lcl.py)
Fields: U10, V10 (m/s), T2 (deg C), Q2 (g/kg), HFX and LH (W/m2, upward positive).
Drift (option --drift, run after the science thread's read of the same-land run A36 is recorded): the same-land member ctrl_p
(F:/CALIFORNIA_wrf_runs/wrf/winter2014_a36b/ctrl_p, 2013-10-01 to 2013-11-08) against ctrl, and pre_w against ctrl, on the days
the three runs share, afternoon only: the per-cell mean difference of each field. The page greys land differences smaller than
the same-land difference in the same cell.
Grid: the page's 5 km grid (valley3d_wind.py), nearest d02 cell within 9 km.
Writes figures/valley_3d/data/wind/wrf_wind.bin (float32; member [ctrl, pre_w], composite, field, row, column), wrf_wind_meta.json;
with --drift also wrf_drift.bin (float32; [ctrl_p - ctrl, pre_w - ctrl], field, row, column) and adds the drift block to the meta.
Prints floor means per member and composite and the pre_w - ctrl difference with a 2,000-sample day resample range.
Run  python src/landscape/valley3d_wrf_wind.py            then, after the A36 read,  python src/landscape/valley3d_wrf_wind.py --drift
"""
import glob
import json
import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd
from netCDF4 import Dataset, chartostring
from rasterio.warp import transform
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import valley3d_wind as W  # noqa: E402  (page 5 km grid and floor cells)

RUN = r"D:\CALIFORNIA_wrf_runs\wrf\winter2014_a28"
DRIFT = r"F:\CALIFORNIA_wrf_runs\wrf\winter2014_a36b\ctrl_p"
OUT = W.OUT
MEMBERS = ["ctrl", "pre_w"]
COMPS = ["afternoon", "predawn", "storm"]
FIELDS = ["U10", "V10", "T2", "Q2", "HFX", "LH"]
T0, T1 = pd.Timestamp("2013-11-01"), pd.Timestamp("2014-03-30 23:00")
NBOOT, SEED = 2000, 20261118


def frames(folder):
    seen = {}
    for f in sorted(glob.glob(os.path.join(folder, "sfc_d02_*"))):
        with Dataset(f) as d:
            for i, s in enumerate(chartostring(d["Times"][:])):
                seen.setdefault(pd.Timestamp(datetime.strptime(str(s), "%Y-%m-%d_%H:%M:%S")), (f, i))
    return seen


def read(idx, times, names):
    byf = {}
    for t in times:
        f, i = idx[t]; byf.setdefault(f, []).append(i)
    out = {v: [] for v in names}
    for f in sorted(byf):
        with Dataset(f) as d:
            for v in names:
                out[v].append(d[v][byf[f]].astype("f8"))
    return {v: np.concatenate(out[v]) for v in names}


def prep(x, v):
    return x - 273.15 if v == "T2" else x * 1000.0 if v == "Q2" else x


# d02 cells: page 5 km grid (nearest within 9 km) and floor cells (d02 centre on a page floor 5 km cell)
wo = sorted(glob.glob(os.path.join(RUN, "ctrl", "wrfout_d02_*")))[0]
with Dataset(wo) as d:
    lat, lon = d["XLAT"][0].astype("f8"), d["XLONG"][0].astype("f8")
xs, ys = transform("EPSG:4326", "EPSG:3310", lon.ravel().tolist(), lat.ravel().tolist())
X, Y = np.array(xs), np.array(ys)
tree = cKDTree(np.c_[X, Y])
GX, GY = np.meshgrid(W.gx, W.gy)
dist, near = tree.query(np.c_[GX.ravel(), GY.ravel()])
gi = np.clip(np.round((X - W.gx[0]) / W.STEP).astype(int), 0, W.WX - 1); gj = np.clip(np.round((W.gy[0] - Y) / W.STEP).astype(int), 0, W.WY - 1)
inside = (X >= W.gx[0]) & (X <= W.gx[-1]) & (Y <= W.gy[0]) & (Y >= W.gy[-1])
floor = inside & W.FLOOR[gj, gi]


def to_page(field2d):
    v = field2d.ravel()[near].astype("f4"); v[dist > 9000] = np.nan
    return v.reshape(W.WY, W.WX)


def main_pair():
    if os.path.exists(os.path.join(OUT, "wrf_wind.bin")):
        raise SystemExit(f"{OUT}/wrf_wind.bin exists; not overwritten")
    idx = {m: frames(os.path.join(RUN, m)) for m in MEMBERS}
    common = sorted(t for t in set(idx["ctrl"]) & set(idx["pre_w"]) if T0 - pd.Timedelta(hours=1) <= t <= T1)
    rr = read(idx["ctrl"], common, ("RAINC", "RAINNC"))
    rain = np.diff(rr["RAINC"] + rr["RAINNC"], axis=0).reshape(len(common) - 1, -1)[:, floor.ravel()].mean(1)
    hours = pd.DatetimeIndex(common[1:])
    sel = {"afternoon": [t for t in hours if t.hour == 21], "predawn": [t for t in hours if t.hour == 12], "storm": list(hours[rain >= 0.1])}
    print("d02 floor cells", int(floor.sum()), "| hours: afternoon", len(sel["afternoon"]), "predawn", len(sel["predawn"]), "storm", len(sel["storm"]))
    A = np.full((2, 3, len(FIELDS), W.WY, W.WX), np.nan, np.float32)
    floor_means, ranges = {}, {}
    rng = np.random.default_rng(SEED)
    for ci, c in enumerate(COMPS):
        S = {m: read(idx[m], sel[c], FIELDS) for m in MEMBERS}
        days = pd.DatetimeIndex(sel[c]).floor("D")
        for fi, v in enumerate(FIELDS):
            fm = {}
            for mi, m in enumerate(MEMBERS):
                x = prep(S[m][v], v)
                A[mi, ci, fi] = to_page(x.mean(0))
                fm[m] = x.reshape(len(x), -1)[:, floor.ravel()].mean(1)   # floor mean per hour
            d = pd.Series(fm["pre_w"] - fm["ctrl"]).groupby(days.values).mean()
            bs = [d.values[rng.integers(0, len(d), len(d))].mean() for _ in range(NBOOT)]
            floor_means[f"{c}_{v}"] = {m: round(float(fm[m].mean()), 3) for m in MEMBERS}
            ranges[f"{c}_{v}"] = [round(float(d.mean()), 3), round(float(np.quantile(bs, .025)), 3), round(float(np.quantile(bs, .975)), 3)]
            print(f"{c:9s} {v:4s} floor ctrl {floor_means[f'{c}_{v}']['ctrl']:9.3f} pre_w {floor_means[f'{c}_{v}']['pre_w']:9.3f} | pre_w - ctrl {ranges[f'{c}_{v}'][0]:+8.3f} [{ranges[f'{c}_{v}'][1]:+.3f}, {ranges[f'{c}_{v}'][2]:+.3f}]")
    A.tofile(os.path.join(OUT, "wrf_wind.bin"))
    meta = {"SOURCE": __doc__.strip(), "order": {"member": MEMBERS, "comp": COMPS, "field": FIELDS}, "grid": {"wx": W.WX, "wy": W.WY},
            "hours": {c: len(sel[c]) for c in COMPS}, "storm_days": len(set(pd.DatetimeIndex(sel["storm"]).floor("D"))),
            "floor_cells_d02": int(floor.sum()), "floor_means": floor_means, "pre_w_minus_ctrl": ranges, "units": {"U10": "m/s", "V10": "m/s", "T2": "deg C", "Q2": "g/kg", "HFX": "W/m2", "LH": "W/m2"},
            "ctrl_isd_dewpoint_depression_bias_k": {"ctrl": 3.43, "pre_w": 0.03}}
    json.dump(meta, open(os.path.join(OUT, "wrf_wind_meta.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print("wrote", OUT)


def main_drift():
    if os.path.exists(os.path.join(OUT, "wrf_drift.bin")):
        raise SystemExit(f"{OUT}/wrf_drift.bin exists; not overwritten")
    idx = {m: frames(os.path.join(RUN, m)) for m in MEMBERS}; idx["ctrl_p"] = frames(DRIFT)
    times = sorted(t for t in set(idx["ctrl"]) & set(idx["pre_w"]) & set(idx["ctrl_p"]) if t.hour == 21)
    print("afternoons shared by ctrl, pre_w and ctrl_p:", len(times), times[0].date() if times else None, "to", times[-1].date() if times else None)
    S = {m: read(idx[m], times, FIELDS) for m in ("ctrl", "pre_w", "ctrl_p")}
    D = np.full((2, len(FIELDS), W.WY, W.WX), np.nan, np.float32); res = {}
    for fi, v in enumerate(FIELDS):
        a, b, c = (prep(S[m][v], v) for m in ("ctrl", "pre_w", "ctrl_p"))
        D[0, fi] = to_page((c - a).mean(0)); D[1, fi] = to_page((b - a).mean(0))
        fl = floor.ravel()
        res[v] = {"same_land_floor": round(float((c - a).reshape(len(times), -1)[:, fl].mean()), 3), "land_floor": round(float((b - a).reshape(len(times), -1)[:, fl].mean()), 3),
                  "same_land_floor_abs_cell": round(float(np.abs((c - a).mean(0)).ravel()[fl].mean()), 3), "land_floor_abs_cell": round(float(np.abs((b - a).mean(0)).ravel()[fl].mean()), 3)}
        print(f"{v:4s} afternoon, floor: same land (ctrl_p - ctrl) {res[v]['same_land_floor']:+.3f}, land (pre_w - ctrl) {res[v]['land_floor']:+.3f}; mean |cell difference| same land {res[v]['same_land_floor_abs_cell']:.3f}, land {res[v]['land_floor_abs_cell']:.3f}")
    D.tofile(os.path.join(OUT, "wrf_drift.bin"))
    mp = os.path.join(OUT, "wrf_wind_meta.json"); meta = json.load(open(mp, encoding="utf-8"))
    meta["drift"] = {"afternoons": len(times), "from": str(times[0].date()), "to": str(times[-1].date()), "order": ["ctrl_p - ctrl", "pre_w - ctrl"], "floor": res}
    json.dump(meta, open(mp, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print("wrote", OUT)


if __name__ == "__main__":
    main_drift() if "--drift" in sys.argv else main_pair()
