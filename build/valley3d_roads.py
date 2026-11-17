"""Major roads for the 3D page (Ali, 2026-11-17: "the road network needs to be shown as well"; "Yes, download it").
Source: US Census TIGER/Line 2023 primary and secondary roads, California (tl_2023_06_prisecroads.zip, 12,741,579 bytes,
www2.census.gov/geo/tiger/TIGER2023/PRISECROADS/), held at Z:/Data/california_central_valley/roads_tiger2023/prisec/.
MTFCC S1100 primary roads (limited-access highways) -> rank 2; S1200 secondary roads (US and state highways) -> rank 1.
Today's network with no build dates: the page shows it as a switch, labelled 2023, not drawn into past years.
Lines clipped to the page grid, simplified to 80 m, in page texture pixels (300 m).
Writes figures/valley_3d/data/roads.json.
WHEN 2026-11-17. STATUS derived (selection and simplification only).
Run  python src/landscape/valley3d_roads.py
"""
import json
import os

import fiona
from rasterio.warp import transform_geom
from shapely.geometry import box, shape

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(REPO, "figures", "valley_3d", "data")
G = json.load(open(os.path.join(DATA, "meta.json"), encoding="utf-8"))["grid"]
SRC = "zip://" + r"Z:\Data\california_central_valley\roads_tiger2023\prisec\tl_2023_06_prisecroads.zip"
grid = box(G["x0"], G["y1"] - G["th"] * G["tex_m"], G["x0"] + G["tw"] * G["tex_m"], G["y1"])
to_px = lambda x, y: [round((x - G["x0"]) / G["tex_m"], 1), round((G["y1"] - y) / G["tex_m"], 1)]
out, n = [], {1: 0, 2: 0}
with fiona.open(SRC) as s:
    for f in s:
        r = 2 if f["properties"]["MTFCC"] == "S1100" else 1
        g = shape(transform_geom(s.crs_wkt, G["crs"], f["geometry"])).intersection(grid).simplify(80)
        parts = [g] if g.geom_type == "LineString" else list(getattr(g, "geoms", []))
        for p in parts:
            if p.geom_type == "LineString" and len(p.coords) > 1:
                out.append({"r": r, "n": f["properties"]["FULLNAME"] or "", "px": [to_px(x, y) for x, y, *_ in p.coords]}); n[r] += 1
json.dump({"SOURCE": __doc__.split("Writes")[0].strip(), "lines": out}, open(os.path.join(DATA, "roads.json"), "w", encoding="utf-8"), separators=(",", ":"))
print("lines", len(out), "| primary", n[2], "| secondary", n[1], "| size", round(os.path.getsize(os.path.join(DATA, "roads.json")) / 1e6, 2), "MB")
