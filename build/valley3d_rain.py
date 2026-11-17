"""Rain for the 3D Central Valley page (Ali, 2026-11-16: "rainfall peaks in months, or season totals, and precipitation
grid on ground, where it falls more").

From PRISM monthly precipitation (4 km, 1895-2025), on a 4 km grid aligned with the page (EPSG:3310, average resampling):
  rain_wy.bin       uint16 little-endian, water-year totals in mm (October to September, labelled by the year it ends),
                    water years 1896-2025, each a rw x rh grid, row 0 north
  rain_meta.json    grid, years, and for every water year the 12 monthly valley-floor means (mm, Oct..Sep) and the
                    floor water-year total; floor = 4 km cells with at least half their area on the page's CVHM2 floor
SOURCE  Z:/Data/california_central_valley/prism_monthly/raw/ppt/prism_ppt_us_25m_<yyyymm>.zip (PRISM Climate Group, OSU)
WHEN    2026-11-16. STATUS derived (resampling and sums; nothing fitted).
Run     python src/landscape/valley3d_rain.py
"""
import glob
import json
import os
import zipfile

import numpy as np
import rasterio
from PIL import Image
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(REPO, "figures", "valley_3d", "data")
OUT = os.path.join(DATA, "rain"); os.makedirs(OUT, exist_ok=True)
G = json.load(open(os.path.join(DATA, "meta.json"), encoding="utf-8"))["grid"]
R_M = 4000.0
RW, RH = int(np.ceil(G["tw"] * G["tex_m"] / R_M)), int(np.ceil(G["th"] * G["tex_m"] / R_M))
T = from_origin(G["x0"], G["y1"], R_M, R_M)
fl = np.array(Image.open(os.path.join(DATA, "floor.png")))
onfloor = ((fl[..., 0].astype(np.int32) + 256 * fl[..., 1].astype(np.int32)) > 0).astype(np.float32)
# floor share per 4 km cell, from the 300 m floor mask resampled
share = np.zeros((RH, RW), np.float32)
reproject(onfloor, share, src_transform=from_origin(G["x0"], G["y1"], G["tex_m"], G["tex_m"]), src_crs=G["crs"], dst_transform=T, dst_crs=G["crs"], resampling=Resampling.average)
floor = share >= 0.5
print(f"rain grid {RW} x {RH} at 4 km | floor cells {int(floor.sum())}")

files = sorted(glob.glob(r"Z:\Data\california_central_valley\prism_monthly\raw\ppt\prism_ppt_us_25m_*.zip"))
month = {}
for path in files:
    ym = os.path.basename(path).split("_")[-1][:6]
    tif = [n for n in zipfile.ZipFile(path).namelist() if n.endswith((".tif", ".bil"))][0]
    with rasterio.open("/vsizip/" + path.replace("\\", "/") + "/" + tif) as src:
        a = np.full((RH, RW), np.nan, np.float32)
        reproject(rasterio.band(src, 1), a, dst_transform=T, dst_crs=G["crs"], resampling=Resampling.average, src_nodata=src.nodata, dst_nodata=np.nan)
    month[ym] = a
    if ym.endswith("12"):
        print("read", ym[:4], flush=True)

years = [y for y in range(1896, 2026) if all(f"{y - (1 if m >= 10 else 0)}{m:02d}" in month for m in list(range(10, 13)) + list(range(1, 10)))]
stack = np.zeros((len(years), RH, RW), np.uint16)
meta = {"SOURCE": __doc__.split("WHEN")[0].strip(), "WHEN": "2026-11-16", "STATUS": "derived",
        "grid": {"crs": G["crs"], "x0": G["x0"], "y1": G["y1"], "cell_m": R_M, "rw": RW, "rh": RH}, "years": years,
        "months": ["Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep"], "floor_cells": int(floor.sum()), "wy": {}}
for n, y in enumerate(years):
    ms = [month[f"{y - (1 if m >= 10 else 0)}{m:02d}"] for m in list(range(10, 13)) + list(range(1, 10))]
    tot = np.nansum(ms, axis=0)
    tot[np.all([np.isnan(a) for a in ms], axis=0)] = np.nan
    stack[n] = np.where(np.isfinite(tot), np.clip(np.round(tot), 1, 65535), 0).astype(np.uint16)
    fm = [round(float(np.nanmean(a[floor])), 1) for a in ms]
    meta["wy"][str(y)] = {"floor_months_mm": fm, "floor_total_mm": round(sum(fm), 1)}
stack.astype("<u2").tofile(os.path.join(OUT, "rain_wy.bin"))
json.dump(meta, open(os.path.join(OUT, "rain_meta.json"), "w", encoding="utf-8"), separators=(",", ":"))
tots = [meta["wy"][str(y)]["floor_total_mm"] for y in years]
print("water years", years[0], "-", years[-1], "| floor total median", float(np.median(tots)), "mm, min", min(tots), "max", max(tots))
print("rain_wy.bin", os.path.getsize(os.path.join(OUT, "rain_wy.bin")) // 1024, "KB")
