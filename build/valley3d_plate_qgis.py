"""QGIS route for the 3D page's print plates (Ali, 2026-11-18: "Both routes").

Reads a request file saved by the page's Plate dialog ("QGIS spec": years, season, extent in EPSG:3310, plate code, title,
subtitle), writes a spec in the grammar of C:/Users/PC/Projects/map-plates/tools/plate.py (the renderer of the 25 plates in
figures/map_plates), renders it headless in QGIS 3.44 and copies the plate to figures/map_plates/page/.

Layers
  land of a year 1874-2019  CVHM2 land use dissolved by class group, the rule of src/landscape/era_plates_prep.py (its GROUPS are
                            imported; that script is not edited or run). Existing lu_<year>.gpkg files are reused; a missing
                            year is written to data/central_valley/page/lu_<year>.gpkg.
  land before 1900          the CVHMP pre-1900 raster used by plate 3 v2 (data/central_valley/ba/rb_hist1900_100m.tif).
  dams                      Global Dam Watch (data/central_valley/dams_cv.gpkg) built by the year, split at 1,000 million m3.
  base                      relief, ocean, floor, states and towns as on the era plates (copied from specs/era_3_1945.json).
Colours: the page's land palette for the request's season (one palette for page and plates, 2026-11-18).
Not drawn: page overlays (wind, wells, rain and the rest); they are drawn by the page's PNG and PDF route. plate.py lays out A3
only, so an A4 request is rendered on A3 and the note says so.
  python src/landscape/valley3d_plate_qgis.py <request.json>
"""
import json
import os
import shutil
import subprocess
import sys

import geopandas as gpd
from shapely.geometry import box

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
from era_plates_prep import CLASS2GROUP, GROUPS, LU  # noqa: E402

MP = r"C:\Users\PC\Projects\map-plates"
CV = os.path.join(MP, "data", "central_valley")
PAGE = os.path.join(CV, "page")
OUT_COPY = os.path.join(REPO, "figures", "map_plates", "page")
QRUN = os.path.join(MP, "tools", "qrun.bat")

# the page's land palette, [winter, summer] (figures/valley_3d/working.html, PAL and CROP first shades)
PAL = {"water": ["#3f6f8a", "#466f80"], "wet_native": ["#46685f", "#5c7a3c"], "dry_native": ["#8eab5c", "#b8a057"],
       "rice": ["#c9bf72", "#d9cf3a"], "orchard": ["#b8a25a", "#c9a227"], "field": ["#d9c86e", "#f2d64b"],
       "pasture": ["#cfc878", "#e6dc6a"], "idle": ["#d8ceaa", "#eadfb4"], "irr_unspec": ["#d9c86e", "#ecd04a"], "urban": ["#5b5d63", "#5b5d63"]}
NAMES = {"water": "Open water", "wet_native": "Wet native cover (marsh, flooded land)", "dry_native": "Dry native cover (grassland, savanna)",
         "rice": "Rice", "orchard": "Orchards, vines, citrus", "field": "Field, row, grain and truck crops", "pasture": "Pasture and hay",
         "idle": "Idle or fallow", "irr_unspec": "Irrigated, crop not classed (CVHM2 class 19)", "urban": "Towns"}
# pre-1900 raster classes (plate 3 v2) in the page's colours: value, page key, label
HIST = [(8, "wet_native", "Wetland (tule marsh)"), (2, "water", "Open water"), (6, "riparian", "Riparian forest"),
        (5, "flood", "Other floodplain habitat"), (7, "oak", "Valley and foothill hardwood"), (4, "dry_native", "Grassland"), (1, "alkali", "Alkali desert scrub")]
PAL_X = {"riparian": ["#465a36", "#3c5a2c"], "flood": ["#5f8a62", "#a9a368"], "oak": ["#7f9450", "#a39a62"], "alkali": ["#b4b394", "#d2c8a6"]}


def land_layer(year):
    """Path (relative to MP) of the dissolved land-use layer of a year and its group areas (km2)."""
    era = os.path.join(CV, "era", f"lu_{year}.gpkg")
    path = era if os.path.exists(era) else os.path.join(PAGE, f"lu_{year}.gpkg")
    if not os.path.exists(path):
        g = gpd.read_file(LU, columns=[f"LU_{year}", "active"]).to_crs(3310)
        g = g[g["active"] == 1].copy()
        g["group"] = g[f"LU_{year}"].fillna(0).astype(int).map(CLASS2GROUP)
        d = g.dropna(subset=["group"])[["group", "geometry"]].dissolve(by="group").reset_index()
        os.makedirs(PAGE, exist_ok=True)
        d.to_file(path, driver="GPKG")
        print("wrote", path)
    d = gpd.read_file(path)
    areas = {r.group: r.geometry.area / 1e6 for r in d.itertuples()}
    return os.path.relpath(path, MP), areas


def dams_layer(year):
    path = os.path.join(PAGE, f"dams_{year}.gpkg")
    if not os.path.exists(path):
        d = gpd.read_file(os.path.join(CV, "dams_cv.gpkg")).to_crs(3310)
        d["storage_mcm"] = d["CAP_MCM"].where(d["CAP_MCM"] > 0)
        d["built"] = d["YEAR_DAM"].where(d["YEAR_DAM"] > 0)
        sel = d[d["built"] <= year]
        os.makedirs(PAGE, exist_ok=True)
        sel.to_file(path, driver="GPKG")
    d = gpd.read_file(path)
    return os.path.relpath(path, MP), len(d), float(d["storage_mcm"].fillna(0).sum()) / 1000


def main(req_path):
    R = json.load(open(req_path, encoding="utf-8"))
    assert R.get("kind") == "valley3d_plate_request", "not a page plate request"
    s = 1 if R["season"] == "summer" else 0
    tmpl = json.load(open(os.path.join(MP, "specs", "era_3_1945.json"), encoding="utf-8"))
    base = {l["id"]: l for l in tmpl["layers"] if l["id"] in ("hs_base", "ocean", "hs_over", "floor", "states", "towns")}
    stem = "page_" + "".join(ch if ch.isalnum() else "_" for ch in f"{R['plate']}_{'_'.join(map(str, R['years']))}").lower().strip("_")
    os.makedirs(PAGE, exist_ok=True)
    ext = os.path.join(PAGE, f"extent_{stem}.gpkg")
    gpd.GeoDataFrame({"id": [1]}, geometry=[box(*R["extent_3310"])], crs=3310).to_file(ext, driver="GPKG")
    layers = [base["hs_base"], base["ocean"], {"id": "view", "name": "View", "path": os.path.relpath(ext, MP), "style": {"kind": "fill", "fill": "#00000000", "stroke": "#00000000", "width": 0}}]
    panels, groups, cols, rows = [], [], [], {}
    for y in R["years"]:
        ids = []
        if y == "pre1900":
            layers.append({"id": "hist", "name": "Before 1900", "type": "raster", "path": "data/central_valley/ba/rb_hist1900_100m.tif",
                           "style": {"kind": "paletted", "opacity": 0.95, "classes": [[v, (PAL.get(k) or PAL_X[k])[s], lab] for v, k, lab in HIST]}})
            ids, label, sub = ["hist"], "Before 1900", "CVHMP pre-1900 vegetation (CSU Chico)"
            groups.append({"name": "Before 1900 (CVHMP)", "layers": ["hist"]})
            cols.append("Before 1900"); rows.setdefault("_note", []).append("Before 1900: drawn from the CVHMP map; areas on the page's PNG plate.")
        else:
            y = int(y); path, areas = land_layer(y)
            for k in GROUPS:
                lid = f"lu_{k}_{y}"
                layers.append({"id": lid, "name": NAMES[k], "path": path, "filter": f"\"group\" = '{k}'", "legend": NAMES[k],
                               "style": {"kind": "fill", "fill": PAL[k][s], "stroke": "#00000000", "width": 0}})
                ids.append(lid)
            dpath, nd, st = dams_layer(y)
            for big in (True, False):
                lid = f"dams_{'big' if big else 'small'}_{y}"
                layers.append({"id": lid, "name": "Dam", "path": dpath, "filter": "\"storage_mcm\" >= 1000" if big else "\"storage_mcm\" < 1000 OR \"storage_mcm\" IS NULL",
                               "legend": "Dam ≥ 1,000 million m³" if big else "Dam, smaller or unknown",
                               "style": {"kind": "point", "shape": "triangle", "fill": "#1f2a30" if big else "#ffffff", "stroke": "#ffffff" if big else "#1f2a30", "size": 3.4 if big else 2.0, "width": 0.4}})
                ids.append(lid)
            label, sub = str(y), f"CVHM2 land use of {y}, dams built by then, {R['season']}"
            if not any(g_["name"].startswith("Land") for g_ in groups):
                groups.append({"name": f"Land (CVHM2), {R['season']} colours", "layers": [f"lu_{k}_{y}" for k in GROUPS]})
                groups.append({"name": "Dams", "layers": [f"dams_big_{y}", f"dams_small_{y}"], "column_break": True})
            cols.append(str(y))
            wet = areas.get("water", 0) + areas.get("wet_native", 0); farm = sum(areas.get(k, 0) for k in ("rice", "orchard", "field", "pasture", "idle", "irr_unspec"))
            for name, v in (("Wet land and open water", wet), ("Dry native land", areas.get("dry_native", 0)), ("Farmed land", farm), ("Towns", areas.get("urban", 0)), ("Dams (storage, km³)", None)):
                rows.setdefault(name, []).append(f"{nd} ({st:.1f})" if v is None else f"{v:,.0f}")
        panels.append({"label": label, "sublabel": sub, "layers": ["hs_base", "ocean"] + ids + ["hs_over", "floor", "states", "towns"]})
    layers += [base["hs_over"], base["floor"], base["states"], base["towns"]]
    tbl_rows = [[k] + v for k, v in rows.items() if k != "_note" and len(v) == len(cols)] or None
    note = (f"Drawn from the 3D page's request ({os.path.basename(req_path)}): the page's 2D view, its season colours. "
            "Page layers such as the wind are drawn on the page's PNG and PDF plates only. CVHM2 is a model input map; class 19 is shown as irrigated, crop not classed (inference). "
            "GDW lacks replaced dams." + (" Requested on A4; this renderer lays out A3." if R.get("paper") == "a4" else "") + " " + " ".join(rows.get("_note", [])))
    title = R["title"] if R.get("overlay", "none") == "none" else f"{R['title']} (land and dams only)"   # page layers are drawn by the PNG and PDF route
    spec = {"plate": R["plate"], "series": "Central Valley in 3D, print plates", "title": title, "title_lines": 2.2 if len(title) < 90 else 3.2,
            "subtitle": R["subtitle"], "subtitle_lines": max(2.2, len(R["subtitle"]) / 70 + 1), "accent": "#0c5560", "crs": "EPSG:3310",
            "crs_label": "NAD83 / California Albers (EPSG:3310), equal-area", "extent_layer": "view", "extent_buffer_km": 0, "grid_deg": 1,
            "scale_km": 50, "scale_segments": 2, "scale_in_map": True, "legend_title": "Legend", "author": "Ali Bin Shahid", "out": f"plates/{stem}",
            "sources_height": 56, "layers": layers, "panels": panels, "legend_groups": groups,
            "sources": ["USGS CVHM2 land use v3, 1874-2019 (cells of about 2.6 km², model active area); CSU Chico CVHMP pre-1900 vegetation; Global Dam Watch v1.0",
                        "Base: USGS WBD HUC12 floor mask; Copernicus GLO-90 relief (© DLR e.V., © Airbus); Natural Earth towns and coast; U.S. Census states"],
            "notes": note.strip()}
    if tbl_rows:
        spec["table"] = {"title": "The floor in numbers, km² (CVHM2 model area)", "header": ["Group"] + cols, "rows": tbl_rows}
    sp = os.path.join(MP, "specs", f"{stem}.json")
    json.dump(spec, open(sp, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("spec", sp)
    r = subprocess.run(["cmd", "/c", QRUN, "-u", os.path.join(MP, "tools", "plate.py"), sp], cwd=os.path.join(MP, "tools"), capture_output=True, text=True)
    print(r.stdout[-1500:], r.stderr[-1500:])
    os.makedirs(OUT_COPY, exist_ok=True)
    for ext_ in (".pdf", ".png", "_preview.png", ".qgz", ".qpt"):
        src = os.path.join(MP, "plates", stem + ext_)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(OUT_COPY, stem + ext_)); print("copied", stem + ext_)
    shutil.copy2(sp, os.path.join(OUT_COPY, stem + ".json"))


if __name__ == "__main__":
    main(sys.argv[1])
