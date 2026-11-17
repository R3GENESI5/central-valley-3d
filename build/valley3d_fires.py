"""Fire years on the 3D Central Valley page grid (300 m, EPSG:3310).

Writes figures/valley_3d/data/fires.png: R, G, B, A = the four most recent fire years at each pixel (year - 1800,
0 = none), R the latest. A pixel burnt by any part of a perimeter of that year (centre inside) takes the year.
SOURCE  Z:/Data/california_central_valley/fire_perimeters/calfire_perimeters_1878_2025_domain_gen50m.geojson
        (valley3d_fire_fetch.py, 20,802 perimeters, simplified 0.0005 deg); grid from figures/valley_3d/data/meta.json.
WHEN    2026-11-16. STATUS derived (rasterising; nothing fitted).
Run     py -3.12 src/landscape/valley3d_fires.py
"""
import collections
import json
import os

import numpy as np
from PIL import Image
from rasterio.features import rasterize
from rasterio.transform import from_origin
from rasterio.warp import transform_geom

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(REPO, "figures", "valley_3d", "data")
G = json.load(open(os.path.join(OUT, "meta.json"), encoding="utf-8"))["grid"]
T = from_origin(G["x0"], G["y1"], G["tex_m"], G["tex_m"])
shape = (G["th"], G["tw"])
src = r"Z:\Data\california_central_valley\fire_perimeters\calfire_perimeters_1878_2025_domain_gen50m.geojson"
by_year = collections.defaultdict(list)
skipped = 0
for f in json.load(open(src, encoding="utf-8"))["features"]:
    y = f["properties"].get("YEAR_")
    if not y or not f.get("geometry"):
        skipped += 1; continue
    by_year[int(y)].append(transform_geom("EPSG:4326", G["crs"], f["geometry"]))
ch = np.zeros((4,) + shape, np.uint8)
area_km2 = {}
for y in sorted(by_year):
    m = rasterize(((g, 1) for g in by_year[y]), out_shape=shape, transform=T, fill=0, dtype="uint8").astype(bool)
    ch[1:][:, m] = ch[:-1][:, m]          # shift older years back
    ch[0][m] = y - 1800
    area_km2[y] = round(float(m.sum()) * G["tex_m"] ** 2 / 1e6)
Image.fromarray(np.moveaxis(ch, 0, -1), "RGBA").save(os.path.join(OUT, "fires.png"))
burnt = ch[0] > 0
print("years", min(by_year), "-", max(by_year), "| perimeters skipped (no year or geometry)", skipped)
print("pixels ever burnt", int(burnt.sum()), "=", round(burnt.sum() * 0.09), "km2 | pixels with more than 4 fires (oldest dropped)", int((ch[3] > 0).sum()))
dec = collections.Counter()
for y, a in area_km2.items():
    dec[y // 10 * 10] += a
print("burnt area on the grid by decade, km2:", dict(sorted(dec.items())))
