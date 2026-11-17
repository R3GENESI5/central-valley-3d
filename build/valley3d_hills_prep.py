"""Change years for the hills and mountains of the 3D Central Valley page (Ali, 2026-11-16: "this still does not look like
there is proper change in the landscape above, i.e. the mountains").

Writes figures/valley_3d/data/hills.png on the page's 300 m grid (EPSG:3310), values as year - 1800 (0 = none):
  R  stand establishment year = 2006 - forest stand age (Pan et al. 2011, North American forest age, 1 km, ca. 2006)
  G  forest loss year (Hansen GFC v1.11 lossyear, 30 m), where at least 25 % of the 300 m pixel lost cover; modal year
  B  first NLCD year (1985-2025 series on disk) in which the pixel's majority class is developed (21-24);
     1985 means developed by 1985 (the page falls back to rough start years before then)
  A  255
SOURCE  Z:/Data/california_central_valley: forest_age_pan2011/ca04_usak06_1km.tif, hansen_gfc/Hansen_GFC-2023-v1.11_lossyear_*,
        nlcd_annual/Annual_NLCD_LndCov_<year>_CU_C1V2.zip; grid from figures/valley_3d/data/meta.json.
WHEN    2026-11-16. STATUS derived (reprojection and per-pixel rules; nothing fitted).
Run     py -3.12 src/landscape/valley3d_hills_prep.py
"""
import glob
import json
import os

import numpy as np
import rasterio
from PIL import Image
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
Z = r"Z:\Data\california_central_valley"
OUT = os.path.join(REPO, "figures", "valley_3d", "data")
G = json.load(open(os.path.join(OUT, "meta.json"), encoding="utf-8"))["grid"]
TW, TH = G["tw"], G["th"]
T = from_origin(G["x0"], G["y1"], G["tex_m"], G["tex_m"])
CRS = G["crs"]

# stand establishment year from forest age (age in years as of about 2006)
age = np.full((TH, TW), np.nan, np.float32)
with rasterio.open(os.path.join(Z, "forest_age_pan2011", "ca04_usak06_1km.tif")) as src:
    reproject(rasterio.band(src, 1), age, dst_transform=T, dst_crs=CRS, resampling=Resampling.nearest, src_nodata=src.nodata, dst_nodata=np.nan)
estab = np.where(np.isfinite(age) & (age > 0), np.clip(2006 - age, 1801, 2006), 0)
print("forest age pixels", int((estab > 0).sum()), "establishment years p5-p95", np.percentile(estab[estab > 0], [5, 50, 95]))

# Hansen loss: fraction of 30 m pixels lost, and the modal loss year among those lost
frac = np.zeros((TH, TW), np.float32); mode = np.zeros((TH, TW), np.float32)
for path in sorted(glob.glob(os.path.join(Z, "hansen_gfc", "Hansen_GFC-2023-v1.11_lossyear_*.tif"))):
    with rasterio.open(path) as src:
        b = src.bounds
        if b.top < 33 or b.bottom > 43 or b.right < -125 or b.left > -116:
            continue
        a = src.read(1)
        f = np.zeros((TH, TW), np.float32); m = np.zeros((TH, TW), np.float32)
        reproject((a > 0).astype(np.float32), f, src_transform=src.transform, src_crs=src.crs, dst_transform=T, dst_crs=CRS, resampling=Resampling.average, src_nodata=None, dst_nodata=0)
        reproject(a.astype(np.float32), m, src_transform=src.transform, src_crs=src.crs, dst_transform=T, dst_crs=CRS, resampling=Resampling.mode, src_nodata=0, dst_nodata=0)
        frac = np.maximum(frac, f); mode = np.where(m > 0, m, mode)
        print("hansen", os.path.basename(path))
loss = np.where((frac >= 0.25) & (mode > 0), 2000 + mode, 0)
print("loss pixels", int((loss > 0).sum()), "by year", np.unique(loss[loss > 0], return_counts=True))

# first developed year from the NLCD series
dev_first = np.zeros((TH, TW), np.int32)
years = sorted(int(os.path.basename(p).split("_")[3]) for p in glob.glob(os.path.join(Z, "nlcd_annual", "Annual_NLCD_LndCov_*_CU_C1V2.zip")))
for y in years:
    path = "/vsizip/" + os.path.join(Z, "nlcd_annual", f"Annual_NLCD_LndCov_{y}_CU_C1V2.zip").replace("\\", "/") + f"/Annual_NLCD_LndCov_{y}_CU_C1V2.tif"
    cls = np.zeros((TH, TW), np.uint8)
    with rasterio.open(path) as src:
        reproject(rasterio.band(src, 1), cls, dst_transform=T, dst_crs=CRS, resampling=Resampling.mode, src_nodata=250, dst_nodata=0)
    dev = (cls >= 21) & (cls <= 24)
    dev_first = np.where((dev_first == 0) & dev, y, dev_first)
    print("nlcd", y, "developed pixels", int(dev.sum()))

enc = lambda v: np.where(v > 0, np.clip(v - 1800, 1, 254), 0).astype(np.uint8)
img = np.dstack([enc(estab), enc(loss), enc(dev_first), np.full((TH, TW), 255, np.uint8)])
Image.fromarray(img, "RGBA").save(os.path.join(OUT, "hills.png"))
print("hills.png", os.path.getsize(os.path.join(OUT, "hills.png")) // 1024, "KB; NLCD years", years)
