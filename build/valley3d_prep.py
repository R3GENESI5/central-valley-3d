"""Data for the 3D Central Valley page (Ali, 2026-11-16: realistic depiction on land cover, before and after with dates;
"Yes, go" to the year-by-year three.js build).

Writes figures/valley_3d/data/:
  dem.png      terrain at 1 km (EPSG:3310): elevation + 500 m = R * 256 + G (m)
  floor.png    300 m grid: R, G = CVHM2 cell index + 1 (little-endian), B = pre-1900 natural type (TYPES), A = 255
  cover.png    300 m grid: R = NLCD 2025 class, G = reservoir: year built - 1800 of the dam a water patch belongs to
               (0 = natural or none), B = 1 on ocean, A = 255
  meta.json    grid, CVHM2 groups and yearly codes (from figures/slider_map/data.json), dams, events, dated canals,
               Tulare 1862 stage, sources
SOURCE  figures/slider_map/data.json (CVHM2 v3 annual land use 1874-2019 grouped as in era_plates_prep.py, GDW dams,
        events of design/LEVEL0_Q2_CHRONOLOGY.md, dated canals, Tulare 1862 stage); Z:/Data/california_central_valley:
        dem_cop90 (Copernicus GLO-90), nlcd_annual 2025, historic_vegetation/pre1900 (CSU Chico CVHMP).
WHEN    2026-11-16. STATUS derived (reprojection, rasterising, lookups; nothing fitted).
Run     py -3.12 src/landscape/valley3d_prep.py
"""
import glob
import io
import json
import os

import fiona
import numpy as np
import rasterio
from PIL import Image
from rasterio.features import rasterize
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject, transform_geom
from scipy import ndimage

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
Z = r"Z:\Data\california_central_valley"
OUT = os.path.join(REPO, "figures", "valley_3d", "data")
os.makedirs(OUT, exist_ok=True)
SLIDER = json.load(open(os.path.join(REPO, "figures", "slider_map", "data.json"), encoding="utf-8"))
F = SLIDER["frame"]; U = F["unit_m"]
MARGIN = 60000.0
X0, Y0 = F["x0"] - MARGIN, F["y0"] - MARGIN
X1, Y1 = F["x0"] + F["w"] * U + MARGIN, F["y0"] + F["h"] * U + MARGIN
TEX, MESH = 300.0, 1000.0
TW, TH = int(round((X1 - X0) / TEX)), int(round((Y1 - Y0) / TEX))
MW, MH = int(round((X1 - X0) / MESH)) + 1, int(round((Y1 - Y0) / MESH)) + 1
T_TEX = from_origin(X0, Y1, TEX, TEX)
CRS = "EPSG:3310"
TYPES = ["none", "grassland", "oak", "riparian", "marsh", "water", "alkali", "floodplain", "chaparral"]
CVHMP_TYPE = {"GRASSLAND": 1, "VALLEY/FOOTHILL HARDWOOD": 2, "RIPARIAN": 3, "WETLAND": 4, "AQUATIC": 5,
              "ALKALI DESERT SCRUB": 6, "OTHER FLOODPLAIN HABITAT": 7, "CHAPARRAL": 8}
to_px = lambda xu, yu: [round((F["x0"] + xu * U - X0) / TEX, 1), round((Y1 - (F["y0"] + yu * U)) / TEX, 1)]
print(f"texture {TW} x {TH} at {TEX:.0f} m, mesh {MW} x {MH} at {MESH:.0f} m")

# ---- terrain: Copernicus GLO-90 tiles averaged onto the 1 km mesh nodes
dem = np.full((MH, MW), np.nan, np.float32)
T_MESH = from_origin(X0 - MESH / 2, Y1 + MESH / 2, MESH, MESH)
for path in sorted(glob.glob(os.path.join(Z, "dem_cop90", "Copernicus_DSM_COG_30_*_DEM.tif"))):
    with rasterio.open(path) as src:
        part = np.full((MH, MW), np.nan, np.float32)
        reproject(rasterio.band(src, 1), part, dst_transform=T_MESH, dst_crs=CRS, resampling=Resampling.average, dst_nodata=np.nan)
        m = np.isfinite(part); dem[m] = part[m]
dem = np.nan_to_num(dem, nan=0.0)
enc = np.clip(np.round(dem + 500), 0, 65535).astype(np.uint16)
Image.fromarray(np.dstack([(enc >> 8).astype(np.uint8), (enc & 255).astype(np.uint8), np.zeros_like(enc, np.uint8)])).save(os.path.join(OUT, "dem.png"))
print("dem", float(dem.min()), float(dem.max()))

# ---- NLCD 2025 at 300 m (majority)
nlcd = np.zeros((TH, TW), np.uint8)
with rasterio.open("/vsizip/" + os.path.join(Z, "nlcd_annual", "Annual_NLCD_LndCov_2025_CU_C1V2.zip").replace("\\", "/") + "/Annual_NLCD_LndCov_2025_CU_C1V2.tif") as src:
    reproject(rasterio.band(src, 1), nlcd, dst_transform=T_TEX, dst_crs=CRS, resampling=Resampling.mode, src_nodata=250, dst_nodata=0)
print("nlcd classes", np.unique(nlcd))

# ---- ocean: zero-elevation nodes outside NLCD land, upsampled to the texture
ocean_mesh = dem <= 0.5
oy = np.clip((np.arange(TH) + 0.5) * TEX / MESH, 0, MH - 1).round().astype(int)
ox = np.clip((np.arange(TW) + 0.5) * TEX / MESH, 0, MW - 1).round().astype(int)
ocean = ocean_mesh[np.ix_(oy, ox)] & ((nlcd == 0) | (nlcd == 11))

# ---- reservoirs: each NLCD water patch off the floor takes the year of the nearest GDW dam within 4 km
cell_idx = np.zeros((TH, TW), np.int32)
C = SLIDER["corners"]
shapes = []
for i in range(len(C) // 8):
    ring = [(F["x0"] + C[i * 8 + 2 * k] * U, F["y0"] + C[i * 8 + 2 * k + 1] * U) for k in range(4)]
    shapes.append(({"type": "Polygon", "coordinates": [ring + [ring[0]]]}, i + 1))
cell_idx = rasterize(shapes, out_shape=(TH, TW), transform=T_TEX, fill=0, dtype="int32")
water = (nlcd == 11) & ~ocean & (cell_idx == 0)
lab, n = ndimage.label(water)
dams = SLIDER["dams"]
dam_px = np.array([to_px(*d["p"]) for d in dams]) if dams else np.zeros((0, 2))
res = np.zeros((TH, TW), np.uint8)
objs = ndimage.find_objects(lab)
dist_px = 4000 / TEX
for k, sl in enumerate(objs, start=1):
    ys, xs = np.nonzero(lab[sl] == k)
    if ys.size < 4:
        continue
    ys = ys + sl[0].start; xs = xs + sl[1].start
    near = None
    for j, (px, py) in enumerate(dam_px):
        if px < xs.min() - dist_px or px > xs.max() + dist_px or py < ys.min() - dist_px or py > ys.max() + dist_px:
            continue
        d = np.min((xs - px) ** 2 + (ys - py) ** 2) ** 0.5
        if d <= dist_px and (near is None or d < near[0]):
            near = (d, dams[j]["built"])
    if near and near[1]:
        res[ys, xs] = int(near[1]) - 1800
print("water patches", n, "reservoir pixels", int((res > 0).sum()))

# ---- pre-1900 natural type on the floor (CVHMP), higher RANK drawn last
recs = []
with fiona.open(os.path.join(Z, "historic_vegetation", "pre1900", "data", "commondata", "pre_1900_final", "pre_1900_final_prj.shp")) as s:
    src_crs = s.crs
    for f in s:
        t = CVHMP_TYPE.get(f["properties"]["CVHMP"])
        if t:
            recs.append((f["properties"]["RANK"] or 0, t, dict(f["geometry"])))
recs.sort(key=lambda r: r[0])
natural = rasterize(((transform_geom(src_crs, CRS, g), t) for _, t, g in recs), out_shape=(TH, TW), transform=T_TEX, fill=0, dtype="uint8")
print("natural type pixels", {TYPES[t]: int((natural == t).sum()) for t in range(1, len(TYPES))})

floor = np.dstack([(cell_idx & 255).astype(np.uint8), (cell_idx >> 8).astype(np.uint8), natural, np.full((TH, TW), 255, np.uint8)])
Image.fromarray(floor, "RGBA").save(os.path.join(OUT, "floor.png"))
cover = np.dstack([nlcd, res, ocean.astype(np.uint8), np.full((TH, TW), 255, np.uint8)])
Image.fromarray(cover, "RGBA").save(os.path.join(OUT, "cover.png"))

# ---- meta
L = SLIDER["layers"]
def lines(feats):
    out = []
    for f in feats:
        for part in f["g"]:
            rings = part if isinstance(part[0][0], list) else [part]
            for r in rings:
                out.append([to_px(*p) for p in r])
    return out
meta = {
    "SOURCE": __doc__.split("WHEN")[0].strip(), "WHEN": "2026-11-16", "STATUS": "derived",
    "grid": {"crs": CRS, "x0": X0, "y1": Y1, "tex_m": TEX, "mesh_m": MESH, "tw": TW, "th": TH, "mw": MW, "mh": MH},
    "types": TYPES, "groups": SLIDER["groups"], "years": SLIDER["years"], "land": SLIDER["land"], "areas": SLIDER["areas"],
    "cells": len(C) // 8,
    "dams": [{"px": to_px(*d["p"]), "name": d["name"], "river": d.get("river"), "built": d["built"], "mcm": d.get("mcm")} for d in dams],
    "events": [{**{k: e[k] for k in ("year", "end", "text", "place", "tag")}, "px": to_px(*e["p"])} for e in SLIDER["events"]],
    "canals": [{"name": f.get("GNIS_Name"), "built": f["built"], "px": lines([f])} for f in L["canals_dated"]],
    "tulare1862": lines(L["tulare1862"]),
    "sources": {**SLIDER["sources"], "terrain": "Copernicus GLO-90 DEM", "cover_hills": "USGS Annual NLCD Collection 1 v2, 2025",
                "pre1900_types": "CSU Chico GIC CVHMP pre-1900 historic vegetation, CVHMP class, rasterised at 300 m by RANK"},
}
json.dump(meta, open(os.path.join(OUT, "meta.json"), "w", encoding="utf-8"), separators=(",", ":"))
print("meta", os.path.getsize(os.path.join(OUT, "meta.json")) // 1024, "KB;", len(meta["events"]), "events,", len(meta["dams"]), "dams,", len(meta["canals"]), "canals")
for f in ("dem.png", "floor.png", "cover.png"):
    print(f, os.path.getsize(os.path.join(OUT, f)) // 1024, "KB")

# ---- cell corners in texture pixels, for drawing fields and roads (added 2026-11-16)
meta["corners_px"] = [round(v, 1) for i in range(len(C) // 8) for k in range(4) for v in to_px(C[i * 8 + 2 * k], C[i * 8 + 2 * k + 1])]
json.dump(meta, open(os.path.join(OUT, "meta.json"), "w", encoding="utf-8"), separators=(",", ":"))
print("meta with corners", os.path.getsize(os.path.join(OUT, "meta.json")) // 1024, "KB")
