"""River and creek network for the 3D page (Ali, 2026-11-17: "the river/stream network is a little unclear").

Sources: NHD flowlines already held for the map plates (nhd_streams.gpkg, EPSG:3310; FType 460 stream/river, 558
artificial path through rivers and lakes; floor and lower foothills), and Natural Earth 10 m rivers of North America
for the mountains outside the NHD extent. Present-day network, drawn in every year of the page.
Ranks: 2 = named river (GNIS name contains "River") or a Natural Earth river; 1 = other named creek or slough;
unnamed flowlines are left out. Lines simplified to 60 m, in page texture pixels (300 m).
Writes figures/valley_3d/data/rivers.json.
WHEN 2026-11-17. STATUS derived (selection and simplification only).
Run  python src/landscape/valley3d_rivers.py
"""
import json
import os

import fiona
from rasterio.warp import transform_geom
from shapely.geometry import box, mapping, shape

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(REPO, "figures", "valley_3d", "data")
G = json.load(open(os.path.join(DATA, "meta.json"), encoding="utf-8"))["grid"]
MP = r"C:\Users\PC\Projects\map-plates\data\central_valley"
grid = box(G["x0"], G["y1"] - G["th"] * G["tex_m"], G["x0"] + G["tw"] * G["tex_m"], G["y1"])
to_px = lambda x, y: [round((x - G["x0"]) / G["tex_m"], 1), round((G["y1"] - y) / G["tex_m"], 1)]

def lines(geom):
    g = geom.intersection(grid).simplify(60)
    parts = [g] if g.geom_type == "LineString" else list(getattr(g, "geoms", []))
    return [[to_px(x, y) for x, y, *_ in p.coords] for p in parts if p.geom_type == "LineString" and len(p.coords) > 1]

out = []
with fiona.open(os.path.join(MP, "nhd_streams.gpkg")) as s:
    nhd_box = box(*s.bounds)
    for f in s:
        nm = f["properties"]["GNIS_Name"]
        if not nm:
            continue
        r = 2 if "River" in nm else 1
        for l in lines(shape(f["geometry"])):
            out.append({"r": r, "n": nm, "px": l})
ne = 0
with fiona.open(os.path.join(MP, "rivers_north_america.gpkg")) as s:
    for f in s:
        g = shape(transform_geom(s.crs_wkt, G["crs"], f["geometry"])).difference(nhd_box)
        if g.is_empty:
            continue
        for l in lines(g):
            out.append({"r": 2, "n": f["properties"].get("name") or "", "px": l}); ne += 1
json.dump({"SOURCE": __doc__.split("Writes")[0].strip(), "lines": out}, open(os.path.join(DATA, "rivers.json"), "w", encoding="utf-8"), separators=(",", ":"))
print("lines", len(out), "| rivers", sum(o["r"] == 2 for o in out), "| creeks", sum(o["r"] == 1 for o in out), "| Natural Earth", ne,
      "| size", round(os.path.getsize(os.path.join(DATA, "rivers.json")) / 1e6, 2), "MB")
