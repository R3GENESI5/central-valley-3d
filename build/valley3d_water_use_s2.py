"""Stage 2 data for the "Where the water goes" support layer (Ali, 2026-11-20: "build stage 2 and stage 3"; plan in
design/WATER_LAYER_STAGES_2_3.md). Writes figures/valley_3d/data/water_use2.json:

  districts  receiving districts drawn as outlines: DWR i03 water-district polygons (texture pixels, simplified to 300 m)
             joined to delivery records by exact name after normalising (case, punctuation, parentheses, "- Ag"/"- M&I",
             WD/ID/WSD/MWC spelled out, "X, City of" -> "city of x"); each with yearly delivered acre-feet per project
             (cvp: USBR CVO tables 22-28 contractor rows 1998-2025; swp: Bulletin 132-25 Table B-5B 1962-2024) and the
             CVP maximum contract quantity (USBR latest-water-contractors.pdf). Unmatched names are listed in the file.
  swp        Table B-5B yearly totals: in-valley (San Joaquin Valley area + Feather River area) and out of the valley
             (South Bay area, North Bay, Central Coast, Southern California), acre-feet. Rows from 2025 are shaded
             (projected) in the bulletin and are not used.
  cities     urban retail suppliers in the Sacramento River, San Joaquin River and Tulare Lake hydrologic regions:
             monthly potable water produced minus water sold to other suppliers (the one production column reported in every year 2014-2026) and yearly demand split (gallons), from SWRCB monthly_combined_dataset.csv, placed
             at the centroid of the supplier's SABL service area(s) (joined on water system ID)
  sewage     wastewater plants inside the grid from the SWRCB volumetric annual reports: yearly effluent by destination
             (acre-feet, reporting years 2019-2025) and, where an eSMR plant has the same name or its outfall lies within 4 km, daily effluent flow (MGD)
SOURCE  Z:/Data/california_central_valley/water_use/* (downloaded 2026-11-20; records raw/tests/2026-11-20_water-use-stage2-downloads)
        and imports/usbr_cvo/raw (CVO tables).
WHEN    2026-11-20. STATUS derived (parsing, joins, unit sums; nothing fitted).
Run     python -P src/landscape/valley3d_water_use_s2.py
"""
import json
import os
import re
import sys

import fitz
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "src", "level0"))
import parse_cvp_deliveries as P  # noqa: E402

W = r"Z:\Data\california_central_valley\water_use"
CVO = r"Z:\Data\california_central_valley\imports\usbr_cvo\raw"
OUT = os.path.join(REPO, "figures", "valley_3d", "data", "water_use2.json")
G = json.load(open(os.path.join(REPO, "figures", "valley_3d", "data", "meta.json"), encoding="utf-8"))["grid"]
GRID = box(G["x0"], G["y1"] - G["th"] * G["tex_m"], G["x0"] + G["tw"] * G["tex_m"], G["y1"])
px = lambda x, y: [round((x - G["x0"]) / G["tex_m"], 1), round((G["y1"] - y) / G["tex_m"], 1)]
log = []

ABBR = [(r"\bwd\b", "water district"), (r"\bid\b", "irrigation district"), (r"\bwsd\b", "water storage district"),
        (r"\bmwc\b", "mutual water company"), (r"\bcwd\b", "county water district"), (r"\bpud\b", "public utility district"),
        (r"\bud\b", "utility district"), (r"\bco\b", "company"), (r"\brd\b", "reclamation district"), (r"\bdist\b", "district"),
        (r"\bwc\b", "water company"), (r"\bwcd\b", "water conservation district"), (r"\bmud\b", "municipal utility district"), (r"#", "")]


def norm(s):
    s = str(s).lower().strip()
    s = re.sub(r"\(.*$", "", s)                      # parentheses, including ones cut off at the line end
    s = re.sub(r"\s+-\s+(ag|m&i|m & i|in|out)\b.*$", "", s)
    m = re.match(r"^(.*),\s*(city|county|town) of$", s.strip())
    if m:
        s = f"{m.group(2)} of {m.group(1)}"
    s = s.replace("&", " and ").replace("-", " ")
    s = re.sub(r"[^a-z0-9 ]", " ", s)
    for a, b in ABBR:
        s = re.sub(a, b, s)
    return re.sub(r"\s+", " ", s).strip()


# ---------------------------------------------------------------- CVP contractor rows, per table per year
def tfile(t, y):
    yy = y % 100
    if y <= 2000:
        return os.path.join(CVO, f"tab{t}{yy:02d}.txt")
    if y <= 2010:
        return os.path.join(CVO, f"tab{t}{yy:02d}.pdf")
    return os.path.join(CVO, f"table_{t}_{y}.pdf")


cvp = {}   # norm name -> {"label", "years": {y: af}}
for t in (22, 23, 25, 26, 27, 28):
    for y in range(1998, 2026):
        fn = tfile(t, y)
        if not os.path.exists(fn):
            continue
        rows = P.rows_txt(fn) if fn.endswith(".txt") else P.rows_pdf(fn)
        block, ok_rows = [], []
        for lab, m, tot in rows:
            l = lab.strip().lower()
            if l.startswith("total"):
                if l == "total" and block:
                    s = sum(r[2] or 0 for r in block)
                    if tot is not None and abs(s - tot) <= 2:
                        ok_rows += block
                    else:
                        log.append(f"cvp t{t} {y}: block of {len(block)} rows sums {s:.0f} vs printed {tot}")
                block = []
                continue
            if tot is None or "pumping" in l or "plus" in l or "wasteway" in l:
                continue
            block.append((lab, m, tot))
        if block:
            log.append(f"cvp t{t} {y}: {len(block)} rows after the last Total not used")
        for lab, m, tot in ok_rows:
            k = norm(lab)
            e = cvp.setdefault(k, {"label": re.sub(r"\s+-\s+(Ag|M&I).*$", "", lab.strip()), "years": {}})
            e["years"][y] = e["years"].get(y, 0) + round(tot)

# CVP maximum contract quantity per contractor (sum over its contracts)
cq = {}
doc = fitz.open(os.path.join(W, "cvp_contracts", "latest-water-contractors.pdf"))
for pg in doc:
    for tb in pg.find_tables().tables:
        for r in tb.extract():
            if not r or not r[0] or r[0].strip().lower().startswith(("water service", "sacramento river settlement", "total")):
                continue
            try:
                q = float(str(r[6]).replace(",", ""))
            except (TypeError, ValueError, IndexError):
                continue
            k = norm(r[0].replace("\n", " "))
            cq[k] = cq.get(k, 0) + q

# ---------------------------------------------------------------- SWP Table B-5B (sheets 1-4 = pages 79-82)
b5 = fitz.open(os.path.join(W, "swp_b132", "B132-25-Appendix-B-032626a.pdf"))
SHEETS = {79: ["Napa", "Solano", "_nb", "Alameda County Flood Control and Water Conservation District Zone 7", "Alameda County Water District",
               "Santa Clara Valley Water District", "_sb", "San Luis Obispo", "Santa Barbara", "_cc"],
          80: ["Dudley Ridge Water District", "Empire West Side Irrigation District", "Kern County Water Agency (M&I)",
               "Kern County Water Agency (Ag)", "Kings", "Oak Flat Water District", "Tulare Lake Basin Water Storage District", "_sjv"],
          81: ["AVEK", "Coachella", "Crestline", "Desert", "Littlerock", "Mojave", "Palmdale", "San Bernardino", "San Gabriel", "San Gorgonio"],
          82: ["Santa Clarita", "Metropolitan", "Ventura", "_sc", "Yuba City", "Butte", "Plumas", "_fr", "_future", "_grand"]}
swp_rows = {}
for pg, cols in SHEETS.items():
    t = [l.strip() for l in b5[pg].get_text().split("\n")]
    for i, l in enumerate(t):
        if re.fullmatch(r"(19[6-9]\d|20[0-2]\d)", l) and int(l) <= 2024:
            vals = t[i + 1:i + 1 + len(cols)]
            try:
                nums = [float(v.replace(",", "")) for v in vals]
            except ValueError:
                continue
            swp_rows.setdefault(int(l), {}).update(dict(zip(cols, nums)))
swp = {}
for y, r in sorted(swp_rows.items()):
    if "_grand" not in r or "_sjv" not in r:
        continue
    sjv_sum = sum(r[c] for c in SHEETS[80][:-1])
    if abs(sjv_sum - r["_sjv"]) > 2:
        log.append(f"swp {y}: SJV contractors {sjv_sum:.0f} vs printed {r['_sjv']:.0f}")
    inv = r["_sjv"] + r.get("_fr", 0)
    swp[y] = {"in_valley_af": round(inv), "out_valley_af": round(r["_grand"] - inv), "south_bay_af": round(r.get("_sb", 0)),
              "southern_ca_af": round(r.get("_sc", 0) + sum(r.get(c, 0) for c in SHEETS[81])), "grand_af": round(r["_grand"]),
              "by": {c: round(r[c]) for c in SHEETS[80][:-1] + ["Yuba City", "Butte", "Plumas"] if c in r}}
SWP_DIST = {"Dudley Ridge Water District": ["Dudley Ridge Water District"], "Empire West Side Irrigation District": ["Empire West Side Irrigation District"],
            "Kern County Water Agency": ["Kern County Water Agency (M&I)", "Kern County Water Agency (Ag)"], "Oak Flat Water District": ["Oak Flat Water District"],
            "Tulare Lake Basin Water Storage District": ["Tulare Lake Basin Water Storage District"]}

# ---------------------------------------------------------------- district outlines
wd = gpd.read_file("zip://" + os.path.join(W, "boundaries", "i03_WaterDistricts_shp.zip")).to_crs(G["crs"])
wd = wd[wd.intersects(GRID)].copy()
wd["k"] = wd["AGENCYNAME"].map(norm)
# a delivery name joins a district when the normalised names are equal, or when the delivery name (at least 12 characters,
# often cut off at the PDF line end) is the start of exactly one district name; the join kind is kept in the file
dkeys = sorted(set(wd["k"]))
link = {}   # district key -> list of (cvp key, how)
for k in cvp:
    if k in dkeys:
        link.setdefault(k, []).append((k, "exact"))
    elif len(k) >= 12:
        hits = [d for d in dkeys if d.startswith(k)]
        if len(hits) == 1:
            link.setdefault(hits[0], []).append((k, "prefix"))
districts, matched = [], set()
for k, g in wd.groupby("k"):
    rec = {}
    if k in link:
        yrs = {}
        for ck, how in link[k]:
            for y, v in cvp[ck]["years"].items():
                yrs[y] = yrs.get(y, 0) + v
            matched.add(ck)
        rec["cvp"] = yrs; rec["join"] = sorted({how for _, how in link[k]}); rec["cvp_names"] = sorted({cvp[ck]["label"] for ck, _ in link[k]})
    for dname, cols in SWP_DIST.items():
        if norm(dname) == k:
            rec["swp"] = {y: sum(swp[y]["by"].get(c, 0) for c in cols) for y in swp}
    if not rec:
        continue
    geom = g.geometry.union_all().intersection(GRID).simplify(300)
    if geom.is_empty:
        continue
    parts = []
    for p in (geom.geoms if geom.geom_type == "MultiPolygon" else [geom]):
        if p.geom_type != "Polygon" or p.area < 1e6:
            continue
        parts.append([[px(x, y) for x, y in p.exterior.coords]] + [[px(x, y) for x, y in r.coords] for r in p.interiors])
    if not parts:
        continue
    c = geom.representative_point()
    districts.append({"name": g["AGENCYNAME"].iloc[0], "area_km2": round(geom.area / 1e6, 1), "c": px(c.x, c.y),
                      "contract_af": round(cq[k]) if k in cq else None, "rings": parts, **rec})
unmatched = sorted(v["label"] for k, v in cvp.items() if k not in matched)

# ---------------------------------------------------------------- cities
u = pd.read_csv(os.path.join(W, "urban_use", "monthly_combined_dataset.csv"), low_memory=False)
u = u[u["HYDRO_REGION"].isin(["Sacramento River", "San Joaquin River", "Tulare Lake"])].copy()
u["m"] = pd.to_datetime(u["REPORT_PERIOD_START_DATE"]).dt.strftime("%Y-%m")
sab = gpd.read_file("zip://" + os.path.join(W, "boundaries", "sabl_public_230322.zip")).to_crs(G["crs"])
sab_by = {str(i).strip(): g for i, g in sab.groupby("SABL_PWSID")}
cities = []
for org, g in u.groupby("ORG_ID"):
    ids = {s.strip() for v in g["WATER_SYSTEM_ID"].dropna() for s in str(v).replace(";", ",").split(",")}
    geoms = [sab_by[i].geometry.union_all() for i in ids if i in sab_by]
    if not geoms:
        log.append(f"city {g['SUPPLIER_NAME'].iloc[0]}: no SABL area for {sorted(ids)[:3]}"); continue
    geo = gpd.GeoSeries(geoms).union_all()
    if not geo.intersects(GRID):
        continue
    c = geo.representative_point()
    mon = g.groupby("m")["POTABLE_SUPPLY_MINUS_SOLD_GAL"].sum(min_count=1).dropna()
    # screen (added 2026-11-20 after Nevada ID reported -1.25e13 gallons in one month): a month is kept only when the figure is
    # zero or more and at most 2,000 gallons per resident per day (about ten times typical use); dropped months are logged
    pop = g.groupby("m")["POP_REPORT_PERIOD"].max().reindex(mon.index)
    bad = (mon < 0) | (mon > 2000 * 31 * pop.fillna(np.inf).where(pop > 0, np.inf))
    for k in mon.index[bad]:
        log.append(f"city {g['SUPPLIER_NAME'].iloc[0]} {k}: {mon[k]:.4g} gallons, population {pop.get(k)}: month not used")
    mon = mon[~bad]
    yr = g.assign(y=g["m"].str[:4].astype(int)).groupby("y")
    split = {int(y): {k: round(float(x[col].sum()), -3) for k, col in (("res", "POTABLE_DEMAND_RES_GAL"), ("cii", "POTABLE_DEMAND_CII_GAL"),
             ("irr", "POTABLE_DEMAND_IRR_GAL"), ("ag", "POTABLE_DEMAND_AG_GAL"))} for y, x in yr}
    cities.append({"name": g["SUPPLIER_NAME"].iloc[0], "c": px(c.x, c.y), "area_km2": round(geo.area / 1e6, 1),
                   "monthly_gal": {k: round(float(v), -3) for k, v in mon.items()}, "split_gal": split})

# ---------------------------------------------------------------- sewage plants
fac = pd.read_csv(os.path.join(W, "wastewater_var", "avr_export_facility.csv"), encoding_errors="replace")
eff = pd.read_csv(os.path.join(W, "wastewater_var", "avr_export_effluent.csv"), encoding_errors="replace", low_memory=False)
fac = fac.sort_values("REPORTING YEAR").groupby("GLOBAL ID").last().reset_index()
fg = gpd.GeoDataFrame(fac, geometry=gpd.points_from_xy(fac["LONGITUDE"], fac["LATITUDE"]), crs="EPSG:4326").to_crs(G["crs"])
fg = fg[fg.within(GRID)]
DEST = {"Discharge to Land": "land", "Inland Surface Waters": "river", "Recycled Water Use": "reuse", "Natural Systems": "river",
        "Enclosed bays, estuaries and coastal lagoons, and ocean waters": "bay", "Underground Injection": "land", "Recycled Water Producer": "reuse"}
eff = eff[eff["GLOBAL ID"].isin(fg["GLOBAL ID"])].copy()
eff["dest"] = eff["DISCHARGED TO"].map(DEST).fillna("other")
eff["af"] = pd.to_numeric(eff["DISCHARGE VOLUME"], errors="coerce")   # VAR volumes are acre-feet (data dictionary)
es = pd.read_parquet(os.path.join(W, "wastewater_esmr", "esmr_flow_region5.parquet"))
es = es[(es["units"] == "MGD") & es["location"].astype(str).str.startswith("EFF") & es["latitude"].notna()]
esp = es.groupby("facility_name").agg(lat=("latitude", "median"), lon=("longitude", "median")).reset_index()
eg = gpd.GeoDataFrame(esp, geometry=gpd.points_from_xy(esp["lon"], esp["lat"]), crs="EPSG:4326").to_crs(G["crs"])
sewage = []
for _, f in fg.iterrows():
    e = eff[eff["GLOBAL ID"] == f["GLOBAL ID"]]
    if e.empty:
        continue
    yrs = {int(y): {d: round(float(v), 1) for d, v in x.groupby("dest")["af"].sum().items()} for y, x in e.groupby("REPORTING YEAR")}
    rec = {"name": f["SITE NAME"], "id": f["GLOBAL ID"], "c": px(f.geometry.x, f.geometry.y), "years_af": yrs}
    same = eg[eg["facility_name"].map(norm) == norm(f["SITE NAME"])]   # same name first, else the nearest outfall within 4 km
    dist = eg.distance(f.geometry)
    nm = same["facility_name"].iloc[0] if len(same) else (eg.loc[dist.idxmin(), "facility_name"] if len(dist) and dist.min() <= 4000 else None)
    if nm and not any(r.get("esmr") == nm for r in sewage):
        d = es[es["facility_name"] == nm].groupby("sampling_date")["result"].sum()
        rec["esmr"] = nm
        rec["daily_mgd"] = {str(k)[:10]: round(float(v), 3) for k, v in d.items()}
    sewage.append(rec)

doc = {"SOURCE": "DWR i03 water districts; USBR CVO tables 22-28 and CVP contractor list; DWR Bulletin 132-25 Table B-5B; SWRCB urban "
                 "supplier monthly data, SABL 2023, volumetric annual reports 2019-2025, CIWQS eSMR (Region 5 flow)",
       "WHEN": "2026-11-20", "STATUS": "derived (parsing, joins, unit sums; nothing fitted)", "AF_M3": 1233.48, "GAL_M3": 0.00378541,
       "districts": districts, "swp": swp, "cities": cities, "sewage": sewage, "unmatched_cvp": unmatched, "parse_log": log}
json.dump(doc, open(OUT, "w", encoding="utf-8"), separators=(",", ":"))

print("size", os.path.getsize(OUT) // 1024, "KB")
print("districts", len(districts), "| with cvp", sum("cvp" in d for d in districts), "| with swp", sum("swp" in d for d in districts),
      "| with contract", sum(d["contract_af"] is not None for d in districts))
print("cvp contractors", len(cvp), "| unmatched", len(unmatched)); print("  unmatched:", "; ".join(unmatched))
print("swp years", min(swp), "-", max(swp), "| 2019", {k: v for k, v in swp[2019].items() if k != "by"})
print("cities", len(cities), "| sewage plants", len(sewage), "| with daily", sum("daily_mgd" in s for s in sewage))
print("parse log", len(log)); [print("  ", l) for l in log[:40]]
for nm in ("Westlands Water District", "Glenn-Colusa Irrigation District", "Kern County Water Agency"):
    d = next((d for d in districts if d["name"] == nm), None)
    if d:
        y = 2019; af = (d.get("cvp", {}).get(y, 0) or 0) + (d.get("swp", {}).get(y, 0) or 0)
        print(f"{nm}: {d['area_km2']} km2, 2019 delivered {af:,} AF = {af * 1233.48 / (d['area_km2'] * 1e6) * 1000:.0f} mm, contract {d['contract_af']}")
