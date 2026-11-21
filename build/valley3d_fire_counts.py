"""Times burnt per pixel on the 3D Central Valley page grid, two equal spans: 1950-1987 and 1988-2025 (38 years each).

Writes figures/valley_3d/data/fire_counts.png: R = fire years 1950-1987 at the pixel, G = fire years 1988-2025,
B = 0, A = 255. A pixel counts a year once if any perimeter of that year covers its centre (same rule as
valley3d_fires.py, which keeps only the four latest years and is not changed).
Spans start in 1950 because the CAL FIRE perimeter record before about 1950 is incomplete
(raw/tests/2026-11-16_valley-3d-fires.md).
SOURCE  Z:/Data/california_central_valley/fire_perimeters/calfire_perimeters_1878_2025_domain_gen50m.geojson
        grid from figures/valley_3d/data/meta.json.
WHEN    2026-11-21. STATUS derived (rasterising and counting; nothing fitted).
Run     py -3.12 src/landscape/valley3d_fire_counts.py
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
SPANS = [(1950, 1987), (1988, 2025)]
src = r"Z:\Data\california_central_valley\fire_perimeters\calfire_perimeters_1878_2025_domain_gen50m.geojson"
by_year = collections.defaultdict(list)
skipped = 0
for f in json.load(open(src, encoding="utf-8"))["features"]:
    y = f["properties"].get("YEAR_")
    if not y or not f.get("geometry"):
        skipped += 1; continue
    by_year[int(y)].append(transform_geom("EPSG:4326", G["crs"], f["geometry"]))
cnt = np.zeros((2,) + shape, np.uint16)
for y in sorted(by_year):
    k = next((i for i, (a, b) in enumerate(SPANS) if a <= y <= b), None)
    if k is None:
        continue
    m = rasterize(((g, 1) for g in by_year[y]), out_shape=shape, transform=T, fill=0, dtype="uint8").astype(bool)
    cnt[k][m] += 1
img = np.zeros(shape + (4,), np.uint8)
img[..., 0] = np.minimum(cnt[0], 255); img[..., 1] = np.minimum(cnt[1], 255); img[..., 3] = 255
Image.fromarray(img, "RGBA").save(os.path.join(OUT, "fire_counts.png"))
px_km2 = G["tex_m"] ** 2 / 1e6
print("perimeters skipped (no year or geometry)", skipped)
for k, (a, b) in enumerate(SPANS):
    c = cnt[k]
    print(f"{a}-{b}: pixels burnt >=1 {int((c >= 1).sum())} ({round((c >= 1).sum() * px_km2):,} km2); "
          f">=2 {int((c >= 2).sum())} ({round((c >= 2).sum() * px_km2):,} km2); >=3 {int((c >= 3).sum())} ({round((c >= 3).sum() * px_km2):,} km2); "
          f"max {int(c.max())}; burnt area counting reburns {round(float(c.sum()) * px_km2):,} km2")
