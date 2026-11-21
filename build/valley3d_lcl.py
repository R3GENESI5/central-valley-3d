"""Cloud-base (LCL) maps, then and now, for the 3D Central Valley page, from the A28 WRF pair (Ali, 2026-11-16: "can we
see LCL then vs now" ... "yes build it, caveated").

Rule: the one of src/wrf/a27_interim_chain.py, applied per cell instead of as a floor mean.
  LCL          125 m per K of T2 - Td2; Td2 from Q2 and PSFC (e = q p / (0.622 + 0.378 q); Bolton 1980 inversion).
  storm hours  hours with CTRL floor-mean rain >= 0.1 mm/h (floor = d02 cells whose centre lies on the page's valley
               floor), the same clock hours in PRE-W (paired timing).
  period       2013-11-01 to 2014-03-30, d02 (12 km), members ctrl (today's land) and pre_w (pre-1900 land).
Writes figures/valley_3d/data/lcl/lcl_<member>.png on the page's 1 km mesh grid (nearest d02 cell; value = LCL m, R*256+G,
0 = no value) and lcl_meta.json (storm hours, floor means, the A28 record's floor result and the CTRL bias for the caveat).
Model world, one winter. CTRL is dry against 12 ISD stations (median model - ISD T - Td +3.43 K; PRE-W +0.03 K):
raw/tests/2026-11-15_a28-pair-read-full-season.
WHEN 2026-11-16. STATUS derived (means over model hours; nothing fitted).
Run  python src/landscape/valley3d_lcl.py
"""
import glob
import json
import os
from datetime import datetime

import numpy as np
import pandas as pd
from netCDF4 import Dataset, chartostring
from PIL import Image
from rasterio.warp import transform
from scipy.spatial import cKDTree

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(REPO, "figures", "valley_3d", "data")
OUT = os.path.join(DATA, "lcl"); os.makedirs(OUT, exist_ok=True)
RUN = r"D:\CALIFORNIA_wrf_runs\wrf\winter2014_a28"
G = json.load(open(os.path.join(DATA, "meta.json"), encoding="utf-8"))["grid"]
fl = np.array(Image.open(os.path.join(DATA, "floor.png")))
cell = fl[..., 0].astype(np.int32) + 256 * fl[..., 1].astype(np.int32)
T0, T1 = pd.Timestamp("2013-11-01"), pd.Timestamp("2014-03-30 23:00")


def frames(m):
    seen = {}
    for f in sorted(glob.glob(os.path.join(RUN, m, "sfc_d02_*"))):
        with Dataset(f) as d:
            for i, s in enumerate(chartostring(d["Times"][:])):
                t = pd.Timestamp(datetime.strptime(str(s), "%Y-%m-%d_%H:%M:%S"))
                seen.setdefault(t, (f, i))
    return seen


def td(q, p):
    e = q * p / (0.622 + 0.378 * q) / 100.0
    ln = np.log(np.maximum(e, 1e-3) / 6.112)
    return 273.15 + 243.5 * ln / (17.67 - ln)


def read(m, times, idx, names):
    byf = {}
    for t in times:
        f, i = idx[t]; byf.setdefault(f, []).append(i)
    out = {v: [] for v in names}
    for f in sorted(byf):
        with Dataset(f) as d:
            for v in names:
                out[v].append(d[v][byf[f]].astype("f8"))
    return {v: np.concatenate(out[v]) for v in names}


# d02 cell centres
wo = sorted(glob.glob(os.path.join(RUN, "ctrl", "wrfout_d02_*")))[0]
with Dataset(wo) as d:
    lat, lon = d["XLAT"][0].astype("f8"), d["XLONG"][0].astype("f8")
xs, ys = transform("EPSG:4326", G["crs"], lon.ravel().tolist(), lat.ravel().tolist())
X, Y = np.array(xs).reshape(lat.shape), np.array(ys).reshape(lat.shape)
px = ((X - G["x0"]) / G["tex_m"]).astype(int); py = ((G["y1"] - Y) / G["tex_m"]).astype(int)
ok = (px >= 0) & (px < G["tw"]) & (py >= 0) & (py < G["th"])
floor = np.zeros(lat.shape, bool); floor[ok] = cell[py[ok], px[ok]] > 0
print("d02", lat.shape, "| floor cells", int(floor.sum()))

idx = {m: frames(m) for m in ("ctrl", "pre_w")}
common = sorted(t for t in set(idx["ctrl"]) & set(idx["pre_w"]) if T0 - pd.Timedelta(hours=1) <= t <= T1)
S = {m: read(m, common, idx[m], ("T2", "Q2", "PSFC", "RAINC", "RAINNC")) for m in idx}
rain = np.diff(S["ctrl"]["RAINC"] + S["ctrl"]["RAINNC"], axis=0)[:, floor].mean(1)
storm = rain >= 0.1
hours = pd.DatetimeIndex(common[1:])
print("hours", len(hours), "| CTRL storm hours", int(storm.sum()), "on", len(set(hours[storm].floor("D"))), "days")

# page mesh nodes -> nearest d02 cell (within 9 km)
tree = cKDTree(np.c_[X.ravel(), Y.ravel()])
mx = G["x0"] + np.arange(G["mw"]) * G["mesh_m"]; my = G["y1"] - np.arange(G["mh"]) * G["mesh_m"]
MX, MY = np.meshgrid(mx, my)
dist, near = tree.query(np.c_[MX.ravel(), MY.ravel()])
meta = {"SOURCE": __doc__.split("WHEN")[0].strip(), "WHEN": "2026-11-16", "STATUS": "derived", "encoding": "LCL m = R * 256 + G; 0 = no value",
        "storm_hours": int(storm.sum()), "storm_days": len(set(hours[storm].floor("D"))), "period": "2013-11-01 to 2014-03-30",
        "a28_record_floor_lcl_storm_m": {"ctrl": 879.255, "pre_w": 773.329, "diff": -105.926, "diff_range": [-150.750, -46.273]},
        "ctrl_isd_dewpoint_depression_bias_k": {"ctrl": 3.43, "pre_w": 0.03}, "members": {}}
for m in ("ctrl", "pre_w"):
    s = {v: S[m][v][1:] for v in ("T2", "Q2", "PSFC")}
    lcl = 125.0 * (s["T2"] - td(s["Q2"], s["PSFC"]))[storm].mean(0)
    v = lcl.ravel()[near]; v[dist > 9000] = np.nan
    grid = v.reshape(MX.shape)
    enc = np.where(np.isfinite(grid), np.clip(np.round(grid), 1, 65535), 0).astype(np.uint16)
    Image.fromarray(np.dstack([(enc >> 8).astype(np.uint8), (enc & 255).astype(np.uint8), np.zeros_like(enc, np.uint8)])).save(os.path.join(OUT, f"lcl_{m}.png"))
    meta["members"][m] = {"floor_mean_m": round(float(lcl[floor].mean()), 1), "domain_p5_p50_p95_m": [round(float(x), 0) for x in np.nanpercentile(lcl, [5, 50, 95])]}
    print(m, "floor mean LCL in CTRL storm hours", meta["members"][m]["floor_mean_m"], "m | domain p5/p50/p95", meta["members"][m]["domain_p5_p50_p95_m"])
json.dump(meta, open(os.path.join(OUT, "lcl_meta.json"), "w", encoding="utf-8"), indent=1)
