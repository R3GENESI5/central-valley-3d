"""Delta tower sites for the 3D page's "where water is returned" layer (main session's four additions, 2026-11-17).
Sites exactly as src/level0/m9_towers_rewetted_air.py lists them (WET rewetted, DRY drained); coordinates from the
AmeriFlux site list held at Z:/Data/california_central_valley/ameriflux/amf_site_list_2026-11-14.json. Numbers on the
page come from raw/tests/2026-11-16_m9-towers-rewetted-air (printed output), not from this file.
Writes figures/valley_3d/data/case/case3.json.
WHEN 2026-11-17. STATUS derived (positions only).
Run  python src/landscape/valley3d_case3.py
"""
import json
import os

from rasterio.warp import transform

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(REPO, "figures", "valley_3d", "data")
G = json.load(open(os.path.join(DATA, "meta.json"), encoding="utf-8"))["grid"]
WET = ["US-Myb", "US-Tw1", "US-Tw4", "US-Tw5", "US-Sne", "US-Twt"]   # as in m9_towers_rewetted_air.py
DRY = ["US-Tw3", "US-Bi1", "US-Snf", "US-Bi2", "US-Tw2"]
sites = {r["SITE_ID"]: r for r in json.load(open(r"Z:\Data\california_central_valley\ameriflux\amf_site_list_2026-11-14.json", encoding="utf-8"))}
out = []
for sid in WET + DRY:
    r = sites[sid]; lat, lon = float(r["GRP_LOCATION"]["LOCATION_LAT"]), float(r["GRP_LOCATION"]["LOCATION_LONG"])
    x, y = transform("EPSG:4326", G["crs"], [lon], [lat])
    out.append([round((x[0] - G["x0"]) / G["tex_m"], 1), round((G["y1"] - y[0]) / G["tex_m"], 1), sid, r["SITE_NAME"], "rewetted" if sid in WET else "drained"])
json.dump({"SOURCE": __doc__.split("Writes")[0].strip(), "towers": out}, open(os.path.join(DATA, "case", "case3.json"), "w", encoding="utf-8"), indent=0)
print(len(out), "towers:", [(o[2], o[4]) for o in out])
