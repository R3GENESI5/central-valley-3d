"""Region outlines for the 3D Central Valley page (Ali, 2026-11-16: "the regions, sacramento, central valley, delta,
san joaquin, tulare etc should have boundaries to show their bounds").

Regions (the project's own units):
  Central Valley floor  union of the 458 valley-floor HUC12 units (data/units/huc12_cv.gpkg, valley_floor = True)
  Sacramento Valley     floor units in HUC4 1802, less the Legal Delta
  San Joaquin Valley    floor units in HUC4 1804, less the Legal Delta
  Tulare Basin          floor units in HUC4 1803
  Delta                 DWR Legal Delta boundary (Z:/Data/california_central_valley/boundaries/legal_delta_boundary_dwr.geojson)
Writes figures/valley_3d/data/regions.json: outer rings in page texture pixels (300 m grid), simplified to 150 m, with
area in km2 and a label point.
WHEN 2026-11-16. STATUS derived (unions and differences of the unit polygons; nothing fitted).
Run  py -3.12 src/landscape/valley3d_regions.py
"""
import json
import os

import fiona
from rasterio.warp import transform_geom
from shapely.geometry import mapping, shape
from shapely.ops import unary_union

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(REPO, "figures", "valley_3d", "data")
G = json.load(open(os.path.join(DATA, "meta.json"), encoding="utf-8"))["grid"]
units = {"1802": [], "1803": [], "1804": []}
with fiona.open(os.path.join(REPO, "data", "units", "huc12_cv.gpkg")) as s:
    assert s.crs.to_epsg() == 3310
    n = 0
    for f in s:
        if f["properties"]["valley_floor"]:
            n += 1
            units.setdefault(f["properties"]["huc4"], []).append(shape(f["geometry"]))
print("floor HUC12 units", n, {k: len(v) for k, v in units.items()})
with fiona.open(r"Z:\Data\california_central_valley\boundaries\legal_delta_boundary_dwr.geojson") as s:
    delta = unary_union([shape(transform_geom(s.crs, "EPSG:3310", dict(f["geometry"]))) for f in s]).buffer(0)
floor = unary_union([g for v in units.values() for g in v]).buffer(0)
regions = [
    ("Central Valley floor", floor),
    ("Sacramento Valley", unary_union(units["1802"]).difference(delta)),
    ("Delta", delta),
    ("San Joaquin Valley", unary_union(units["1804"]).difference(delta)),
    ("Tulare Basin", unary_union(units["1803"])),
]
px = lambda x, y: [round((x - G["x0"]) / G["tex_m"], 1), round((G["y1"] - y) / G["tex_m"], 1)]
out = []
for name, geom in regions:
    g = geom.simplify(150)
    parts = [g] if g.geom_type == "Polygon" else list(g.geoms)
    rings = [[px(x, y) for x, y in p.exterior.coords] for p in parts if p.area > 2e6]
    lp = geom.representative_point()
    out.append({"name": name, "area_km2": round(geom.area / 1e6), "rings": rings, "label": px(lp.x, lp.y)})
    print(f"{name}: {geom.area / 1e6:,.0f} km2, {len(rings)} ring(s)")
json.dump({"WHEN": "2026-11-16", "SOURCE": __doc__.split("WHEN")[0].strip(), "regions": out},
          open(os.path.join(DATA, "regions.json"), "w", encoding="utf-8"), separators=(",", ":"))
