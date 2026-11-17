"""Four case layers for the 3D Central Valley page (Ali, 2026-11-16: "Top four").

Writes figures/valley_3d/data/case/:
  case.json       wells then (WSP 1905-14, 4,648 located entries) and now (DWR spring 2015-2024 medians, 5,558 wells) as
                  [px, py, class]; rain-yield stations (A7, 29 with a dated drop) as [px, py, name, year, change];
                  rim rivers (14) with decade natural flow, flow below the dam and share; floor evaporation by season
                  then and now (km3, mm) from the before-and-after plate statistics
  evap_<when>_<season>.png   evaporation in mm per season on the page's 1 km mesh grid (R*256+G, 0 = no value),
                  when = then | now, season = nov_feb | mar_may | jun_sep
SOURCE  C:/Users/PC/Projects/map-plates/data/central_valley/ba/: p2_wsp_wells.gpkg, p2_dwr_wells_spring_2015_2024.gpkg,
        p2_stats.json, p4_{then,now}_<season>_mm_1km.tif, p4_stats.json (plates ba_02 and ba_04, 2026-11-06,
        raw/tests/2026-11-06_before-after-maps-eleven-pairs); figures/slider_map/data.json (stations from
        raw/tests/2026-10-23_a7-event-transients-full-86yr; flow from src/layer4/wb_model.py pairing).
WHEN 2026-11-16. STATUS derived (reprojection and lookups of held results; nothing fitted).
Run  python src/landscape/valley3d_case.py
"""
import json
import os

import fiona
import numpy as np
import rasterio
from PIL import Image
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(REPO, "figures", "valley_3d", "data")
OUT = os.path.join(DATA, "case"); os.makedirs(OUT, exist_ok=True)
BA = r"C:\Users\PC\Projects\map-plates\data\central_valley\ba"
G = json.load(open(os.path.join(DATA, "meta.json"), encoding="utf-8"))["grid"]
S = json.load(open(os.path.join(REPO, "figures", "slider_map", "data.json"), encoding="utf-8"))
F = S["frame"]; U = F["unit_m"]
px = lambda x, y: [round((x - G["x0"]) / G["tex_m"], 1), round((G["y1"] - y) / G["tex_m"], 1)]
pu = lambda xu, yu: px(F["x0"] + xu * U, F["y0"] + yu * U)
CLS = {"under 3 m": 0, "3 to 7.6 m": 1, "over 7.6 m": 2, "flowing": 3, "flowing or at surface": 3}


def wells(path):
    out = []
    with fiona.open(path) as s:
        assert s.crs.to_epsg() == 3310
        for f in s:
            p = f["properties"]
            if p.get("on_floor") is False or p.get("cls") not in CLS:
                continue
            x, y = f["geometry"]["coordinates"][:2]
            out.append(px(x, y) + [CLS[p["cls"]]])
    return out


then = wells(os.path.join(BA, "p2_wsp_wells.gpkg"))
now = wells(os.path.join(BA, "p2_dwr_wells_spring_2015_2024.gpkg"))
count = lambda w: [sum(1 for r in w if r[2] == k) for k in range(4)]
p2 = json.load(open(os.path.join(BA, "p2_stats.json"), encoding="utf-8"))
p4 = json.load(open(os.path.join(BA, "p4_stats.json"), encoding="utf-8"))
case = {
    "SOURCE": __doc__.split("WHEN")[0].strip(), "WHEN": "2026-11-16", "STATUS": "derived",
    "well_classes": ["within 3 m", "3 to 7.6 m", "deeper than 7.6 m", "flowing or at the surface"],
    "wells_then": then, "wells_now": now, "wells_then_counts": count(then), "wells_now_counts": count(now),
    "wells_median_m": {"then_nonflowing": p2.get("wsp_median_depth_m_nonflowing"), "now_on_floor": p2.get("dwr_median_depth_m_on_floor")},
    "stations": [pu(*s["p"]) + [s["name"], s["year"], s["change"]] for s in S["stations"]],
    "rivers": [{"river": r["river"], "site": r["site"], "px": pu(*r["p"]), "decades": r["decades"]} for r in S["flow"]],
    "evap": {s: {k: p4[s][k] for k in ("then_floor_km3", "now_floor_km3", "then_floor_mean_mm", "now_floor_mean_mm")} for s in ("nov_feb", "mar_may", "jun_sep")},
    "evap_now_product": p4["now_product"],
}
T = from_origin(G["x0"] - G["mesh_m"] / 2, G["y1"] + G["mesh_m"] / 2, G["mesh_m"], G["mesh_m"])
for when in ("then", "now"):
    for s in ("nov_feb", "mar_may", "jun_sep"):
        with rasterio.open(os.path.join(BA, f"p4_{when}_{s}_mm_1km.tif")) as src:
            a = np.full((G["mh"], G["mw"]), np.nan, np.float32)
            reproject(rasterio.band(src, 1), a, dst_transform=T, dst_crs=G["crs"], resampling=Resampling.nearest, src_nodata=src.nodata, dst_nodata=np.nan)
        enc = np.where(np.isfinite(a) & (a >= 0), np.clip(np.round(a), 1, 65535), 0).astype(np.uint16)
        Image.fromarray(np.dstack([(enc >> 8).astype(np.uint8), (enc & 255).astype(np.uint8), np.zeros_like(enc, np.uint8)])).save(os.path.join(OUT, f"evap_{when}_{s}.png"))
        print(when, s, "cells", int((enc > 0).sum()), "mean mm", round(float(np.nanmean(a[enc > 0])), 1))
json.dump(case, open(os.path.join(OUT, "case.json"), "w", encoding="utf-8"), separators=(",", ":"))
print("wells then", len(then), case["wells_then_counts"], "| now", len(now), case["wells_now_counts"], "| median m", case["wells_median_m"])
print("stations", len(case["stations"]), "| rivers", len(case["rivers"]), "| evap km3", {s: (v["then_floor_km3"], v["now_floor_km3"]) for s, v in case["evap"].items()})
