"""Peak surface temperature for the 3D Central Valley page (Ali, 2026-11-16: "valley-wide hottest days, go").

Rule, fixed before any number was looked at:
  - Source: MODIS Aqua MYD11A1 v6.1 daytime land surface temperature (overpass about 13:30 local), 1 km, AppEEARS cubes.
  - Good pixel: LST_Day_1km > 0, mandatory QA (QC_Day bits 0-1) = 0 or 1, LST error (bits 6-7) <= 1 (2 K or less).
  - Valley floor: MODIS pixels whose centre falls in a CVHM2 cell of the page (floor.png).
  - Candidate day: 1 June to 30 September, at least 70 % of floor pixels good.
  - The year's 10 valley-wide hottest days: the 10 candidate days with the highest floor-mean LST.
  - Peak map: per pixel, the mean LST over those 10 days where at least 4 of them are good; otherwise no value.
  - Years: 2003-2025 (Aqua's first full summer is 2003).
Writes figures/valley_3d/data/lst/peak_<year>.png on the page's 1 km mesh grid (EPSG:3310; value = round((LST C - 10) * 4),
0 = no value) and figures/valley_3d/data/lst/lst_meta.json (the 10 days per year with floor mean and coverage).
SOURCE  Z:/Data/california_central_valley/modis_lst/ca_cv_myd11a1_*/MYD11A1.061_1km_aid0001.nc
WHEN    2026-11-16. STATUS derived (screening, means, reprojection; nothing fitted).
Run     py -3.12 src/landscape/valley3d_lst_peaks.py
"""
import datetime as dt
import glob
import json
import os

import h5py
import numpy as np
from PIL import Image
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject, transform

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(REPO, "figures", "valley_3d", "data")
OUT = os.path.join(DATA, "lst"); os.makedirs(OUT, exist_ok=True)
G = json.load(open(os.path.join(DATA, "meta.json"), encoding="utf-8"))["grid"]
floor_png = np.array(Image.open(os.path.join(DATA, "floor.png")))
cell = floor_png[..., 0].astype(np.int32) + 256 * floor_png[..., 1].astype(np.int32)
EPOCH = dt.datetime(2000, 1, 1)
COVER_MIN, NDAYS, MINGOOD = 0.70, 10, 4

floor_mask = None
days = {}   # year -> list of (date, floor mean C, coverage, file, index)
files = sorted(glob.glob(r"Z:\Data\california_central_valley\modis_lst\ca_cv_myd11a1_*\MYD11A1.061_1km_aid0001.nc"))
for path in files:
    with h5py.File(path, "r") as f:
        lat, lon = f["lat"][:], f["lon"][:]
        if floor_mask is None:
            LO, LA = np.meshgrid(lon, lat)
            xs, ys = transform("EPSG:4326", G["crs"], LO.ravel(), LA.ravel())
            px = ((np.array(xs) - G["x0"]) / G["tex_m"]).astype(int); py_ = ((G["y1"] - np.array(ys)) / G["tex_m"]).astype(int)
            ok = (px >= 0) & (px < G["tw"]) & (py_ >= 0) & (py_ < G["th"])
            fm = np.zeros(px.size, bool); fm[ok] = cell[py_[ok], px[ok]] > 0
            floor_mask = fm.reshape(LO.shape)
            src_T = from_origin(lon[0] - (lon[1] - lon[0]) / 2, lat[0] - (lat[1] - lat[0]) / 2, lon[1] - lon[0], lat[0] - lat[1])
            print("MODIS grid", lat.size, "x", lon.size, "| floor pixels", int(floor_mask.sum()), flush=True)
        t = f["time"][:]
        for k, sec in enumerate(t):
            d = EPOCH + dt.timedelta(seconds=int(sec))
            if d.year < 2003 or not (6 <= d.month <= 9):
                continue
            lst = f["LST_Day_1km"][k]; qc = f["QC_Day"][k]
            good = (lst > 0) & ((qc & 3) <= 1) & (((qc >> 6) & 3) <= 1)
            gf = good & floor_mask
            cov = gf.sum() / floor_mask.sum()
            if cov >= COVER_MIN:
                days.setdefault(d.year, []).append((d.date().isoformat(), float(lst[gf].mean() * 0.02 - 273.15), round(float(cov), 3), path, k))
    print("scanned", os.path.basename(os.path.dirname(path)), flush=True)

meta = {"SOURCE": __doc__.split("WHEN")[0].strip(), "WHEN": "2026-11-16", "STATUS": "derived", "encoding": "LST C = 10 + value / 4; 0 = no value", "years": {}}
dst_T = from_origin(G["x0"] - G["mesh_m"] / 2, G["y1"] + G["mesh_m"] / 2, G["mesh_m"], G["mesh_m"])
for y in sorted(days):
    top = sorted(days[y], key=lambda r: -r[1])[:NDAYS]
    acc = None; n = None
    for date, fmean, cov, path, k in top:
        with h5py.File(path, "r") as f:
            lst = f["LST_Day_1km"][k].astype(np.float64); qc = f["QC_Day"][k]
        good = (lst > 0) & ((qc & 3) <= 1) & (((qc >> 6) & 3) <= 1)
        v = np.where(good, lst * 0.02 - 273.15, 0.0)
        acc = v if acc is None else acc + v; n = good.astype(int) if n is None else n + good
    peak = np.where(n >= MINGOOD, acc / np.maximum(n, 1), np.nan).astype(np.float32)
    out = np.full((G["mh"], G["mw"]), np.nan, np.float32)
    reproject(peak, out, src_transform=src_T, src_crs="EPSG:4326", dst_transform=dst_T, dst_crs=G["crs"], resampling=Resampling.bilinear, src_nodata=np.nan, dst_nodata=np.nan)
    enc = np.where(np.isfinite(out), np.clip(np.round((out - 10) * 4), 1, 255), 0).astype(np.uint8)
    Image.fromarray(enc, "L").save(os.path.join(OUT, f"peak_{y}.png"))
    fl = peak[floor_mask & np.isfinite(peak)]
    meta["years"][str(y)] = {"days": [{"date": d, "floor_mean_c": round(m, 2), "floor_cover": c} for d, m, c, _, _ in top],
                             "candidate_days": len(days[y]), "floor_peak_mean_c": round(float(fl.mean()), 2) if fl.size else None}
    print(y, "candidates", len(days[y]), "| hottest days", [d for d, *_ in top][:3], "... | floor peak mean", meta["years"][str(y)]["floor_peak_mean_c"], "C", flush=True)
json.dump(meta, open(os.path.join(OUT, "lst_meta.json"), "w", encoding="utf-8"), indent=1)
print("wrote", len(meta["years"]), "years")
