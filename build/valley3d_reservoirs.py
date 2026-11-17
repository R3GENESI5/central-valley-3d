"""Reservoir years for the 3D Central Valley page from GDW reservoir outlines (replaces the nearest-dam rule of
valley3d_prep.py, which left reservoirs whose dam lay more than 4 km from the water, such as Trinity Lake, as water in 1880).

Rewrites the G channel of figures/valley_3d/data/cover.png: for each NLCD 2025 water pixel (class 11) touched by a GDW
reservoir outline, year built - 1800; 0 for natural water. Outlines flagged LAKE_CTRL = "Yes" (a natural lake regulated
by a dam, e.g. Lake Tahoe) are natural and always water.
SOURCE  Z:/Data/california_central_valley/gdw_dams/GDW_v1_0_shp.zip (GDW_reservoirs_v1_0: YEAR_DAM, LAKE_CTRL).
WHEN    2026-11-16. STATUS derived (rasterising; nothing fitted).
Run     py -3.12 src/landscape/valley3d_reservoirs.py
"""
import json
import os

import fiona
import numpy as np
from PIL import Image
from rasterio.features import rasterize
from rasterio.transform import from_origin
from rasterio.warp import transform_geom

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(REPO, "figures", "valley_3d", "data")
G = json.load(open(os.path.join(OUT, "meta.json"), encoding="utf-8"))["grid"]
T = from_origin(G["x0"], G["y1"], G["tex_m"], G["tex_m"])
cover = np.array(Image.open(os.path.join(OUT, "cover.png")))
shapes, natural = [], []
src = "zip://Z:/Data/california_central_valley/gdw_dams/GDW_v1_0_shp.zip!GDW_v1_0_shp/GDW_reservoirs_v1_0.shp"
with fiona.open(src) as s:
    for f in s.filter(bbox=(-124.5, 33.5, -117.0, 42.0)):
        p = f["properties"]; g = transform_geom(s.crs, G["crs"], dict(f["geometry"]))
        if p.get("LAKE_CTRL") == "Yes":
            natural.append(p.get("RES_NAME"))
        elif p.get("YEAR_DAM") and p["YEAR_DAM"] > 1800:
            shapes.append((g, int(p["YEAR_DAM"]) - 1800))
shapes.sort(key=lambda r: r[1])   # where outlines overlap, the later dam is drawn last
yr = rasterize(shapes, out_shape=cover.shape[:2], transform=T, fill=0, all_touched=True, dtype="uint8")
water = cover[..., 0] == 11
before = int((cover[..., 1] > 0).sum())
cover[..., 1] = np.where(water, yr, 0)
Image.fromarray(cover, "RGBA").save(os.path.join(OUT, "cover.png"))
print("reservoir outlines", len(shapes), "| natural lakes kept", natural)
print("reservoir pixels", before, "->", int((cover[..., 1] > 0).sum()), "| water pixels left natural", int((water & (cover[..., 1] == 0)).sum()))
