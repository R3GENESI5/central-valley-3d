"""Water rights for the "Where the water goes" support layer, version b (Ali, 2026-11-20: depiction plan v2 in
design/WATER_LAYER_STAGES_2_3.md; screening changes proposed by the science session the same day; "show me first before adding
it to the main"). Successor of valley3d_water_rights_s3.py (recorded, not edited). Writes figures/valley_3d/data/water_rights_b.json.

Rules, written before any total was read:
  rights    active surface-water POD inside the page grid; U.S. Bureau of Reclamation and DWR rights left out (their water is in
            stages 1-2); groundwater recordations left out.
  types     per right, year and month: DIRECT, STORAGE and "Combined (Direct + Storage)" kept apart (combined is what older
            reports give; it cannot be split). USE rows not used. Years to 2020 from the cleaned demand-analysis set, 2021-2024
            from the raw reported-use archive (cache derived_raw_2021plus_grid_rows.parquet; one row per right, month, type).
  limit     (a) the right's own limit (PODs' rated rate x 365 days + POD storage; for a face value, face value + maximum storage):
            a year above it is "above" and not drawn. (b) a right with no limit on record is checked against its own other
            years: a year above 10 x the median of its reported non-zero years is "above"; a right with fewer than 3 reported
            years and no limit is "unchecked" and not drawn. Drawn, above and unchecked volumes are kept per year.
            (c) added after the first run: rights with no limit on record that pass (b) are drawn but kept outside the totals
            ("nolimit" flag), because statements that repeat one unit error every year pass (b) (S019257 Steamboat Slough
            1,776 million m3 and S003103 White Spring 1,480 million m3 in 2019).
  flags     drawn but outside the totals: power (any Power use; hydropower water returns to the river), owner is a federal
            contractor (CVP tables 22-28) or a state contractor (Bulletin 132 Table B-5B valley and Feather contractors),
            amount repeated by the same owner, source and year (one statement filed twice).
  huc8      HUC 8 watershed of the right's first POD (demand-analysis POD file), for the aggregated circles at the overview zoom.
Units: million m3 per right and year (1 acre-foot = 1233.48 m3).
SOURCE  Z:/Data/california_central_valley/water_use/water_rights/ (records raw/tests/2026-11-20_water-rights-stage3-downloads,
        raw/tests/2026-11-20_water-rights-demand-water-use-download).
WHEN    2026-11-20. STATUS derived (filters, de-duplication, unit sums; nothing fitted).
Run     python -P src/landscape/valley3d_water_rights_s3b.py
"""
import json
import os
import re

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
R = r"Z:\Data\california_central_valley\water_use\water_rights"
OUT = os.path.join(REPO, "figures", "valley_3d", "data", "water_rights_b.json")
G = json.load(open(os.path.join(REPO, "figures", "valley_3d", "data", "meta.json"), encoding="utf-8"))["grid"]
GRID = box(G["x0"], G["y1"] - G["th"] * G["tex_m"], G["x0"] + G["tw"] * G["tex_m"], G["y1"])
px = lambda x, y: [round((x - G["x0"]) / G["tex_m"], 1), round((G["y1"] - y) / G["tex_m"], 1)]
MM3 = 1233.48e-6          # million m3 per acre-foot
CFS_AFY = 723.97
RATE_AFY = {"Cubic Feet per Second": CFS_AFY, "Gallons per Day": 365 / 325851.0, "Gallons per Minute": 1440 * 365 / 325851.0,
            "Acre-feet per Year": 1.0, "Acre-feet": 1.0, "Gallons": 1 / 325851.0}
ABBR = [(r"\bwd\b", "water district"), (r"\bid\b", "irrigation district"), (r"\bwsd\b", "water storage district"),
        (r"\bmwc\b", "mutual water company"), (r"\bcwd\b", "county water district"), (r"\bpud\b", "public utility district"),
        (r"\bud\b", "utility district"), (r"\bco\b", "company"), (r"\brd\b", "reclamation district"), (r"\bdist\b", "district"),
        (r"\bwc\b", "water company"), (r"\bwcd\b", "water conservation district"), (r"\bmud\b", "municipal utility district"), (r"#", ""),
        (r"\bcalif\b", "california")]


def norm(s):
    s = str(s).lower().strip()
    s = re.sub(r"\(.*$", "", s)
    s = re.sub(r"\s+-\s+(ag|m&i|m & i|in|out)\b.*$", "", s)
    m_ = re.match(r"^(.*),\s*(city|county|town) of$", s.strip())
    if m_:
        s = f"{m_.group(2)} of {m_.group(1)}"
    s = s.replace("&", " and ").replace("-", " ")
    s = re.sub(r"[^a-z0-9 ]", " ", s)
    for a, b in ABBR:
        s = re.sub(a, b, s)
    return re.sub(r"\s+", " ", s).strip()


# rights and their PODs inside the grid
pod = pd.read_csv(os.path.join(R, "pod", "ewrims_flat_file_pod-flat-file.csv"), low_memory=False, encoding_errors="replace")
pod = pod[pod["LATITUDE"].notna() & pod["LONGITUDE"].notna()]
pg = gpd.GeoDataFrame(pod, geometry=gpd.points_from_xy(pod["LONGITUDE"], pod["LATITUDE"]), crs="EPSG:4326").to_crs(G["crs"])
pg = pg[pg.within(GRID) & (pg["POD_STATUS"].astype(str) == "Active") & ~pg["POD_TYPE"].astype(str).str.contains("Well|Groundwater", case=False)]
rt = pd.read_csv(os.path.join(R, "demand_analysis", "ewrims_flat_file.csv"), low_memory=False, encoding_errors="replace").drop_duplicates("APPLICATION_NUMBER").set_index("APPLICATION_NUMBER")
appl = sorted(set(pg["APPL_ID"].astype(str)) & set(rt.index))
own = rt.loc[appl, "PRIMARY_OWNER_NAME"].astype(str)
project = own.str.upper().str.contains("BUREAU OF RECLAMATION|DEPARTMENT OF WATER RESOURCES|DEPT OF WATER RESOURCES|U S BUREAU|US BUREAU")
gw = rt.loc[appl, "WATER_RIGHT_TYPE"].astype(str).str.contains("Groundwater", case=False)
keep = sorted(own[~project & ~gw].index)
us = pd.read_csv(os.path.join(R, "demand_analysis", "ewrims_flat_file_use_season.csv"), low_memory=False, encoding_errors="replace")
power = set(us.loc[us["USE_CODE"].astype(str) == "Power", "APPLICATION_NUMBER"].astype(str))
dpod = pd.read_csv(os.path.join(R, "demand_analysis", "ewrims_flat_file_pod.csv"), usecols=["APPL_ID", "HUC_8_NUMBER", "HUC_8_NAME"], low_memory=False, encoding_errors="replace").dropna().drop_duplicates("APPL_ID").set_index("APPL_ID")

# own limit (a)
pk = pod[pod["APPL_ID"].astype(str).isin(keep)].copy()
pk["works_afy"] = pd.to_numeric(pk["DIRECT_DIVERSION_RATE"], errors="coerce").fillna(0) * pk["DIRECT_DIV_RATE_UNIT"].map(RATE_AFY).fillna(0) \
    + pd.to_numeric(pk["STORAGE_AMOUNT"], errors="coerce").fillna(0)
works = pk.drop_duplicates(["APPL_ID", "POD_ID"]).groupby(pk["APPL_ID"].astype(str))["works_afy"].sum()
limit = {}
for a in keep:
    w = float(works.get(a, 0) or 0)
    r = rt.loc[a]; fv = pd.to_numeric(r["FACE_VALUE_AMOUNT"], errors="coerce"); ms = pd.to_numeric(r["MAX_STORAGE"], errors="coerce")
    fa = float(fv) * (CFS_AFY if "Cubic" in str(r["FACE_VALUE_UNITS"]) else 1.0) + (float(ms) if ms and ms > 0 else 0.0) if fv and fv > 0 else None
    lim = [x for x in (w if w > 0 else None, fa) if x]
    limit[a] = min(lim) if lim else None

# monthly amounts by type
TYPES = {"DIRECT": "d", "STORAGE": "s", "Combined (Direct + Storage)": "c"}
cl = pd.read_csv(os.path.join(R, "demand_analysis", "water_use_report.csv"), usecols=["APPL_ID", "YEAR", "MONTH", "AMOUNT", "DIVERSION_TYPE"])
cl = cl[cl["APPL_ID"].isin(set(keep)) & (cl["YEAR"] <= 2020) & cl["DIVERSION_TYPE"].isin(TYPES)]
cl["src"] = "clean"
raw = pd.read_parquet(os.path.join(R, "use_reported", "derived_raw_2021plus_grid_rows.parquet"))
raw = raw[raw["APPL_ID"].isin(set(keep)) & raw["DIVERSION_TYPE"].isin(TYPES)].copy()
raw["AMOUNT"] = pd.to_numeric(raw["AMOUNT"], errors="coerce"); raw["src"] = "raw"
m = pd.concat([cl, raw]).drop_duplicates(["APPL_ID", "YEAR", "MONTH", "DIVERSION_TYPE"])
m["t"] = m["DIVERSION_TYPE"].map(TYPES)
yt = m.pivot_table(index=["APPL_ID", "YEAR"], columns="t", values="AMOUNT", aggfunc="sum").fillna(0)
for k in "dsc":
    if k not in yt:
        yt[k] = 0.0
yt["tot"] = yt["d"] + yt["s"] + yt["c"]
yt["src"] = m.groupby(["APPL_ID", "YEAR"])["src"].first()
yt = yt.reset_index()

# screen (a) and (b)
med = yt[yt["tot"] > 0].groupby("APPL_ID")["tot"].agg(["median", "size"])
def flag(row):
    lim = limit.get(row["APPL_ID"])
    if lim is not None:
        return "ok" if row["tot"] <= lim else "above"
    if row["APPL_ID"] not in med.index or med.loc[row["APPL_ID"], "size"] < 3:
        return "unchecked"
    return "ok" if row["tot"] <= 10 * med.loc[row["APPL_ID"], "median"] else "above"
yt["flag"] = yt.apply(flag, axis=1)
yt = yt[yt["tot"] > 0]
scr = yt.groupby(["YEAR", "flag"])["tot"].agg(["size", "sum"])
screened = {}
for (y, f), v in scr.iterrows():
    screened.setdefault(int(y), {})[f] = {"n": int(v["size"]), "mm3": round(float(v["sum"]) * MM3, 2)}

# owner overlaps: federal (CVP tables) and state (Table B-5B valley and Feather contractors)
W2 = json.load(open(os.path.join(REPO, "figures", "valley_3d", "data", "water_use2.json"), encoding="utf-8"))
cvpn = {k for k in ({norm(n) for d in W2["districts"] for n in d.get("cvp_names", [])} | {norm(n) for n in W2["unmatched_cvp"]}) if len(k) >= 6}
SWPN = [norm(x) for x in ("Dudley Ridge Water District", "Empire West Side Irrigation District", "Kern County Water Agency", "Oak Flat Water District",
                          "Tulare Lake Basin Water Storage District", "City of Yuba City", "County of Butte", "County of Plumas", "Kings County")]
def contractor(owner):
    k = norm(owner)
    if k in cvpn or any(len(c) >= 12 and k.startswith(c) for c in cvpn):
        return "cvp"
    if any(k == c or (len(c) >= 12 and k.startswith(c)) for c in SWPN):
        return "swp"
    return None

ok = yt[yt["flag"] == "ok"]
rights, seen = [], {}
for a, g in ok.groupby("APPL_ID"):
    r = rt.loc[a]; p = pg[pg["APPL_ID"].astype(str) == a]; q = p.geometry.iloc[0]
    years = {}
    for _, x in g.sort_values("YEAR").iterrows():
        years[int(x["YEAR"])] = {"d": round(x["d"] * MM3, 4), "s": round(x["s"] * MM3, 4), "c": round(x["c"] * MM3, 4), "src": x["src"]}
    rec = {"appl": a, "type": str(r["WATER_RIGHT_TYPE"]), "owner": str(r["PRIMARY_OWNER_NAME"]),
           "priority": None if pd.isna(r["PRIORITY_DATE"]) else str(r["PRIORITY_DATE"])[:10],
           "face": None if pd.isna(r["FACE_VALUE_AMOUNT"]) or not r["FACE_VALUE_AMOUNT"] else round(float(r["FACE_VALUE_AMOUNT"]) * (CFS_AFY if "Cubic" in str(r["FACE_VALUE_UNITS"]) else 1.0) * MM3, 4),
           "source": str(r["SOURCE_NAME"]), "c": px(q.x, q.y), "n_pod": int(len(p)),
           "huc8": None if a not in dpod.index else str(dpod.loc[a, "HUC_8_NUMBER"]).split(".")[0], "huc8_name": None if a not in dpod.index else str(dpod.loc[a, "HUC_8_NAME"]),
           "power": a in power, "contractor": contractor(r["PRIMARY_OWNER_NAME"]), "nolimit": limit.get(a) is None, "years": years}
    rights.append(rec)
for rec in sorted(rights, key=lambda x: x["appl"]):
    for y, v in rec["years"].items():
        tot = v["d"] + v["s"] + v["c"]
        key = (norm(rec["owner"]), rec["source"].lower(), y, round(tot, 4))
        if key in seen:
            v["dup_of"] = seen[key]
        else:
            seen[key] = rec["appl"]

# HUC8 groups for the overview: position = mean of its rights' first PODs
hucs = {}
for rec in rights:
    if rec["huc8"]:
        h = hucs.setdefault(rec["huc8"], {"name": rec["huc8_name"], "xy": []}); h["xy"].append(rec["c"])
huc8 = {k: {"name": v["name"], "c": [round(float(np.mean([p[0] for p in v["xy"]])), 1), round(float(np.mean([p[1] for p in v["xy"]])), 1)], "n": len(v["xy"])} for k, v in hucs.items()}

doc = {"SOURCE": "SWRCB eWRIMS PODs, demand-analysis rights, uses and cleaned water use (to 2020); raw reported use 2021-2024",
       "WHEN": "2026-11-20", "STATUS": "derived (filters, de-duplication, unit sums; nothing fitted)", "UNIT": "million m3",
       "rights": rights, "huc8": huc8, "screened": screened}
json.dump(doc, open(OUT, "w", encoding="utf-8"), separators=(",", ":"))

print("size", os.path.getsize(OUT) // 1024, "KB | rights kept", len(keep), "| drawn", len(rights), "| HUC8 groups", len(huc8))
print("flags: power", sum(r["power"] for r in rights), "| cvp owner", sum(r["contractor"] == "cvp" for r in rights), "| swp owner", sum(r["contractor"] == "swp" for r in rights))
for y in (1995, 2000, 2010, 2015, 2019, 2021, 2023):
    inn = [r for r in rights if y in r["years"] and not r["power"] and not r["contractor"] and not r["nolimit"] and not r["years"][y].get("dup_of")]
    nl = sum(sum(r["years"][y][k] for k in "dsc") for r in rights if y in r["years"] and r["nolimit"] and not r["power"] and not r["contractor"]) / 1e3
    pw = sum(sum(r["years"][y][k] for k in "dsc") for r in rights if y in r["years"] and r["power"]) / 1e3
    d = sum(r["years"][y]["d"] for r in inn) / 1e3; s = sum(r["years"][y]["s"] for r in inn) / 1e3; c = sum(r["years"][y]["c"] for r in inn) / 1e3
    out = sum(sum(r["years"][y][k] for k in "dsc") for r in rights if y in r["years"] and (r["power"] or r["contractor"] or r["nolimit"] or r["years"][y].get("dup_of"))) / 1e3
    sc = screened.get(y, {})
    print(f"{y}: in totals direct {d:.2f} | storage {s:.2f} | combined {c:.2f} km3 ({len(inn)} rights) | flagged {out:.2f} (power {pw:.2f}, no limit {nl:.2f}) | above {sc.get('above', {}).get('mm3', 0) / 1e3:.2f} ({sc.get('above', {}).get('n', 0)}) | unchecked {sc.get('unchecked', {}).get('mm3', 0) / 1e3:.2f} ({sc.get('unchecked', {}).get('n', 0)})")
