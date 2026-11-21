"""Dam walls and reservoir sheets for the 3D Central Valley page (Ali, 2026-11-20: "do a dam pass, shouldn't it have a
reservoir behind the wall"; plan agreed "Yes, go").

For every dam in figures/valley_3d/data/meta.json, finds its Global Dam Watch barrier (same build year, nearest point
within 3 km; else the same name), and writes figures/valley_3d/data/dams.json:
  dams: [{i (index in meta.dams), gdw (GDW_ID), hgt_m, len_m, area_km2, elev_masl, use, res: [[ring, ...], ...]}]
  res  is the dam's GDW reservoir outline in texture pixels (300 m grid of the page), one list of rings per polygon part,
       first ring the shore, the rest islands; simplified to 60 m; parts under 0.05 km2 dropped. Empty when the
       reservoir has no outline or the outline is flagged LAKE_CTRL = "Yes" (a natural lake regulated by a dam).
SOURCE  Z:/Data/california_central_valley/gdw_dams/GDW_v1_0_shp.zip (GDW_barriers_v1_0, GDW_reservoirs_v1_0).
WHEN    2026-11-20. STATUS derived (matching, reprojection, simplification; nothing fitted).
Run     py -3.12 src/landscape/valley3d_dam_shapes.py
"""
import json
import os

import fiona
import numpy as np
from rasterio.warp import transform_geom
from shapely.geometry import shape, mapping

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(REPO, "figures", "valley_3d", "data")
META = json.load(open(os.path.join(OUT, "meta.json"), encoding="utf-8"))
G = META["grid"]
ZIP = "zip://Z:/Data/california_central_valley/gdw_dams/GDW_v1_0_shp.zip!GDW_v1_0_shp/"
BBOX = (-124.5, 33.5, -117.0, 42.0)


def px(x, y):
    return ((x - G["x0"]) / G["tex_m"], (G["y1"] - y) / G["tex_m"])


bar = []
with fiona.open(ZIP + "GDW_barriers_v1_0.shp") as s:
    for f in s.filter(bbox=BBOX):
        p = f["properties"]; g = transform_geom(s.crs, G["crs"], dict(f["geometry"]))
        x, y = g["coordinates"][:2]
        bar.append({"id": p["GDW_ID"], "name": (p["DAM_NAME"] or "").lower(), "year": p["YEAR_DAM"], "px": px(x, y),
                    "hgt": p["DAM_HGT_M"], "len": p["DAM_LEN_M"], "area": p["AREA_SKM"], "elev": p["ELEV_MASL"],
                    "use": p["MAIN_USE"], "lakectrl": p["LAKE_CTRL"]})
res = {}
with fiona.open(ZIP + "GDW_reservoirs_v1_0.shp") as s:
    for f in s.filter(bbox=BBOX):
        p = f["properties"]
        if p.get("LAKE_CTRL") == "Yes":
            continue
        res[p["GDW_ID"]] = shape(transform_geom(s.crs, G["crs"], dict(f["geometry"])))

out, unmatched, by = [], [], {"year+3km": 0, "name": 0}
bxy = np.array([b["px"] for b in bar])
for i, d in enumerate(META["dams"]):
    dist = np.hypot(bxy[:, 0] - d["px"][0], bxy[:, 1] - d["px"][1]) * G["tex_m"]
    cand = [k for k in np.argsort(dist) if dist[k] <= 3000 and bar[k]["year"] == d["built"]]
    if cand:
        b = bar[cand[0]]; by["year+3km"] += 1
    else:
        nm = [k for k in np.argsort(dist) if bar[k]["name"] == (d["name"] or "").lower()]
        if not nm:
            unmatched.append(d["name"]); continue
        b = bar[nm[0]]; by["name"] += 1
    rings = []
    poly = res.get(b["id"])
    if poly is not None:
        poly = poly.simplify(60)
        for part in (poly.geoms if poly.geom_type == "MultiPolygon" else [poly]):
            if part.area < 50000:
                continue
            rr = [part.exterior] + list(part.interiors)
            rings.append([[[round(v, 1) for v in px(x, y)] for x, y in r.coords] for r in rr])
    out.append({"i": i, "gdw": b["id"], "hgt_m": b["hgt"], "len_m": b["len"], "area_km2": b["area"],
                "elev_masl": b["elev"], "use": b["use"], "res": rings})

doc = {"SOURCE": "GDW v1.0 barriers and reservoirs (Z:/Data/california_central_valley/gdw_dams/GDW_v1_0_shp.zip)",
       "WHEN": "2026-11-20", "STATUS": "derived (matching, reprojection, simplification)", "dams": out}
json.dump(doc, open(os.path.join(OUT, "dams.json"), "w", encoding="utf-8"), separators=(",", ":"))
print("dams in meta", len(META["dams"]), "| matched", len(out), by, "| unmatched", len(unmatched), unmatched[:10])
print("with reservoir outline", sum(1 for o in out if o["res"]), "| without", sum(1 for o in out if not o["res"]))
print("height m: min", min(o["hgt_m"] or 0 for o in out), "max", max(o["hgt_m"] or 0 for o in out),
      "| missing height", sum(1 for o in out if not o["hgt_m"]), "| missing length", sum(1 for o in out if not o["len_m"]))
print("size", os.path.getsize(os.path.join(OUT, "dams.json")) // 1024, "KB")
for nm in ("Shasta", "Oroville", "Folsom"):
    o = next(o for o in out if META["dams"][o["i"]]["name"] == nm)
    print(nm, {k: o[k] for k in ("gdw", "hgt_m", "len_m", "area_km2", "elev_masl")}, "parts", len(o["res"]),
          "points", sum(len(r) for p in o["res"] for r in p))
