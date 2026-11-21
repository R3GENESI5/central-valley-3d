"""Stage 3 data for the "Where the water goes" support layer: water taken from the rivers under water rights (Ali, 2026-11-20:
"build stage 2 and stage 3"; plan in design/WATER_LAYER_STAGES_2_3.md). Writes figures/valley_3d/data/water_rights.json:

  rights  one entry per water right with at least one active surface-water point of diversion inside the page grid:
          application number, right type, owner, priority date, face value, source, the first POD inside the grid (texture
          pixels) and the count of PODs; yearly diverted acre-feet:
            diverted = Combined (Direct + Storage) where reported that month, else DIRECT + STORAGE (USE rows not added:
            they report water used, which overlaps what was diverted)
          years 1906-2020 from the State Board's cleaned demand-analysis set (water_use_report.csv; its 2021 is incomplete); 2021-2024 from the raw
          eWRIMS reported-use archive (the 7z), one row kept per right, year, month and diversion type (the archive repeats
          the amount once per joined POD/use/party row). Each year carries its source ("clean" or "raw").
  Screen: rights with a Power use are left out (hydropower water returns to the river). A year is drawn only when the
          diverted amount is at most the right's own limit: rated rate of its PODs for 365 days plus POD storage, and for rights
          with a face value no more than face value plus maximum storage. Reports above the limit or with no limit on record
          are counted per year in "screened" and not drawn.
  Left out: rights owned by the U.S. Bureau of Reclamation or the California Department of Water Resources (their water is
  drawn in stages 1 and 2), groundwater recordations, and rights with no diversion reported in any year.
SOURCE  Z:/Data/california_central_valley/water_use/water_rights/ (records raw/tests/2026-11-20_water-rights-stage3-downloads,
        raw/tests/2026-11-20_water-rights-demand-water-use-download).
WHEN    2026-11-20. STATUS derived (filters, de-duplication, unit sums; nothing fitted).
Run     python -P src/landscape/valley3d_water_rights_s3.py
"""
import json
import os
import subprocess

import geopandas as gpd
import pandas as pd
from shapely.geometry import box

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
R = r"Z:\Data\california_central_valley\water_use\water_rights"
OUT = os.path.join(REPO, "figures", "valley_3d", "data", "water_rights.json")
G = json.load(open(os.path.join(REPO, "figures", "valley_3d", "data", "meta.json"), encoding="utf-8"))["grid"]
GRID = box(G["x0"], G["y1"] - G["th"] * G["tex_m"], G["x0"] + G["tw"] * G["tex_m"], G["y1"])
px = lambda x, y: [round((x - G["x0"]) / G["tex_m"], 1), round((G["y1"] - y) / G["tex_m"], 1)]
SEVENZ = r"C:\Program Files\7-Zip\7z.exe"

# points of diversion inside the grid
pod = pd.read_csv(os.path.join(R, "pod", "ewrims_flat_file_pod-flat-file.csv"), low_memory=False, encoding_errors="replace")
pod = pod[pod["LATITUDE"].notna() & pod["LONGITUDE"].notna()]
pg = gpd.GeoDataFrame(pod, geometry=gpd.points_from_xy(pod["LONGITUDE"], pod["LATITUDE"]), crs="EPSG:4326").to_crs(G["crs"])
pg = pg[pg.within(GRID)]
print("PODs inside grid", len(pg), "| status", pg["POD_STATUS"].value_counts().to_dict())
pg = pg[(pg["POD_STATUS"].astype(str) == "Active") & ~pg["POD_TYPE"].astype(str).str.contains("Well|Groundwater", case=False)]

rt = pd.read_csv(os.path.join(R, "demand_analysis", "ewrims_flat_file.csv"), low_memory=False, encoding_errors="replace").drop_duplicates("APPLICATION_NUMBER")
rt = rt.set_index("APPLICATION_NUMBER")
appl = sorted(set(pg["APPL_ID"].astype(str)) & set(rt.index))
own = rt.loc[appl, "PRIMARY_OWNER_NAME"].astype(str).str.upper()
project = own.str.contains("BUREAU OF RECLAMATION|DEPARTMENT OF WATER RESOURCES|DEPT OF WATER RESOURCES|U S BUREAU|US BUREAU")
gw = rt.loc[appl, "WATER_RIGHT_TYPE"].astype(str).str.contains("Groundwater", case=False)
us = pd.read_csv(os.path.join(R, "demand_analysis", "ewrims_flat_file_use_season.csv"), low_memory=False, encoding_errors="replace")
power = set(us.loc[us["USE_CODE"].astype(str) == "Power", "APPLICATION_NUMBER"].astype(str))   # hydropower: the water returns to the river
pw = pd.Series([a in power for a in appl], index=appl)
keep = set(own[~project & ~gw & ~pw].index)
print("rights with an active surface POD in grid", len(appl), "| project rights left out", int(project.sum()), "| groundwater", int(gw.sum()),
      "| power use left out", int((pw & ~project & ~gw).sum()), "| kept", len(keep))

# the right's own limit, acre-feet a year: what its diversion works can carry (rated rate at each POD for 365 days, plus POD
# storage), and for licences and permits no more than the face value (plus maximum storage). A year's report above the limit,
# or of a right with no limit on record, is not drawn and is counted separately.
CFS_AFY = 723.97                                    # 1 cfs for 365 days = 723.97 acre-feet
RATE_AFY = {"Cubic Feet per Second": CFS_AFY, "Gallons per Day": 365 / 325851.0, "Gallons per Minute": 1440 * 365 / 325851.0,
            "Acre-feet per Year": 1.0, "Acre-feet": 1.0, "Gallons": 1 / 325851.0}
pk = pod[pod["APPL_ID"].astype(str).isin(keep)].copy()
pk["works_afy"] = pd.to_numeric(pk["DIRECT_DIVERSION_RATE"], errors="coerce").fillna(0) * pk["DIRECT_DIV_RATE_UNIT"].map(RATE_AFY).fillna(0) \
    + pd.to_numeric(pk["STORAGE_AMOUNT"], errors="coerce").fillna(0)
works = pk.drop_duplicates(["APPL_ID", "POD_ID"]).groupby(pk["APPL_ID"].astype(str))["works_afy"].sum()
limit = {}
for a in keep:
    w = float(works.get(a, 0) or 0)
    r = rt.loc[a]; fv = pd.to_numeric(r["FACE_VALUE_AMOUNT"], errors="coerce"); ms = pd.to_numeric(r["MAX_STORAGE"], errors="coerce")
    fa = None
    if fv and fv > 0:
        fa = float(fv) * (CFS_AFY if "Cubic" in str(r["FACE_VALUE_UNITS"]) else 1.0) + (float(ms) if ms and ms > 0 else 0.0)
    lim = [x for x in (w if w > 0 else None, fa) if x]
    limit[a] = min(lim) if lim else None


def monthly_diverted(df):
    """per right, year, month: Combined where present, else DIRECT + STORAGE"""
    df = df[df["DIVERSION_TYPE"].isin(["DIRECT", "STORAGE", "Combined (Direct + Storage)"])]
    p = df.pivot_table(index=["APPL_ID", "YEAR", "MONTH"], columns="DIVERSION_TYPE", values="AMOUNT", aggfunc="first")
    comb = p.get("Combined (Direct + Storage)")
    ds = p.get("DIRECT", 0).fillna(0) + p.get("STORAGE", 0).fillna(0) if ("DIRECT" in p or "STORAGE" in p) else 0
    out = comb.where(comb.notna(), ds) if comb is not None else ds
    return out.rename("af").reset_index()


cl = pd.read_csv(os.path.join(R, "demand_analysis", "water_use_report.csv"), usecols=["APPL_ID", "YEAR", "MONTH", "AMOUNT", "DIVERSION_TYPE"])
cl = cl[cl["APPL_ID"].isin(keep) & (cl["YEAR"] <= 2020)]   # the cleaned set holds 2021 for 14,716 rights against 30,697 in 2020: 2021 is taken from the raw archive
mc = monthly_diverted(cl)
mc["src"] = "clean"

cols = ["APPL_ID", "YEAR", "MONTH", "AMOUNT", "DIVERSION_TYPE"]
CACHE = os.path.join(R, "use_reported", "derived_raw_2021plus_grid_rows.parquet")   # rows of the 7z for 2021+ and grid rights (no change to values)
if os.path.exists(CACHE):
    raw = pd.read_parquet(CACHE); nrow = "cached"
else:
    proc = subprocess.Popen([SEVENZ, "e", "-so", os.path.join(R, "use_reported", "water-rights-water-use-reported.7z")], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    parts, nrow = [], 0
    allgrid = set(appl)
    for ch in pd.read_csv(proc.stdout, usecols=cols, chunksize=1_000_000, low_memory=False, encoding="utf-8-sig", encoding_errors="replace"):
        nrow += len(ch)
        ch = ch[(ch["YEAR"] >= 2021) & ch["APPL_ID"].isin(allgrid)]
        if len(ch):
            parts.append(ch.drop_duplicates(["APPL_ID", "YEAR", "MONTH", "DIVERSION_TYPE"]))
    proc.wait()
    raw = pd.concat(parts).drop_duplicates(["APPL_ID", "YEAR", "MONTH", "DIVERSION_TYPE"]) if parts else pd.DataFrame(columns=cols)
    raw.to_parquet(CACHE, index=False)
raw = raw[raw["APPL_ID"].isin(keep)].copy()
raw["AMOUNT"] = pd.to_numeric(raw["AMOUNT"], errors="coerce")
mr = monthly_diverted(raw)
mr["src"] = "raw"
print("raw archive rows read", nrow, "| raw 2022+ rows kept", len(raw), "| years", sorted(raw["YEAR"].unique().tolist()))

m = pd.concat([mc, mr])
yr = m.groupby(["APPL_ID", "YEAR"]).agg(af=("af", "sum"), months=("af", "size"), src=("src", "first")).reset_index()
yr["limit"] = yr["APPL_ID"].map(limit)
yr["flag"] = "ok"
yr.loc[yr["limit"].isna(), "flag"] = "nolimit"
yr.loc[yr["limit"].notna() & (yr["af"] > yr["limit"]), "flag"] = "above"
screen = yr[yr["af"] > 0].groupby(["YEAR", "flag"])["af"].agg(["size", "sum"]).unstack("flag").fillna(0)
screened = {int(y): {f: {"n": int(screen.loc[y, ("size", f)]), "af": round(float(screen.loc[y, ("sum", f)]), 1)} for f in ("above", "nolimit") if ("size", f) in screen.columns}
            for y in screen.index}
ok = yr[yr["flag"] == "ok"]
rights = []
for a, g in ok.groupby("APPL_ID"):
    if g["af"].sum() <= 0:
        continue
    r = rt.loc[a]; p = pg[pg["APPL_ID"].astype(str) == a]; q = p.geometry.iloc[0]
    rights.append({"appl": a, "type": str(r["WATER_RIGHT_TYPE"]), "owner": str(r["PRIMARY_OWNER_NAME"]),
                   "priority": None if pd.isna(r["PRIORITY_DATE"]) else str(r["PRIORITY_DATE"])[:10],
                   "face_af": None if pd.isna(r["FACE_VALUE_AMOUNT"]) else float(r["FACE_VALUE_AMOUNT"]), "face_unit": str(r["FACE_VALUE_UNITS"]),
                   "source": str(r["SOURCE_NAME"]), "c": px(q.x, q.y), "n_pod": int(len(p)), "limit_afy": round(limit[a], 1),
                   "years": {int(y): {"af": round(float(v), 1), "src": s} for y, v, s in zip(g["YEAR"], g["af"], g["src"])}})

doc = {"SOURCE": "SWRCB eWRIMS points of diversion; demand-analysis rights and cleaned water use (1906-2020); raw reported use 2021-2024",
       "WHEN": "2026-11-20", "STATUS": "derived (filters, de-duplication, unit sums; nothing fitted)", "rights": rights,
       "screened": screened}
json.dump(doc, open(OUT, "w", encoding="utf-8"), separators=(",", ":"))
print("size", os.path.getsize(OUT) // 1024, "KB | rights drawn", len(rights))
tot = pd.Series({y: sum(r["years"].get(y, {}).get("af", 0) for r in rights) for y in range(2008, 2025)})
print("diverted, all kept rights, km3:", {int(y): round(v * 1233.48e-9, 2) for y, v in tot.items()})
print("by type 2019 km3:", pd.Series({r["type"]: 0 for r in rights}).add(pd.Series([r["years"].get(2019, {}).get("af", 0) * 1233.48e-9 for r in rights], index=[r["type"] for r in rights]).groupby(level=0).sum(), fill_value=0).round(2).to_dict())
top = sorted(rights, key=lambda r: -r["years"].get(2019, {}).get("af", 0))[:10]
for r in top:
    print(f"  {r['appl']} {r['type'][:12]:12} {r['owner'][:40]:40} {r['source'][:28]:28} 2019 {r['years'].get(2019, {}).get('af', 0):,.0f} AF")
print("screened (not drawn) 2019:", screened.get(2019), "| 2023:", screened.get(2023))
print("limit known for", sum(v is not None for v in limit.values()), "of", len(limit), "kept rights")

# ---- overlaps (added 2026-11-20 after the first totals): a right is flagged, drawn but left out of the yearly total, when
#  (a) its owner is a federal (CVP) contractor in USBR CVO tables 22-28 (settlement and exchange contractors report the same
#      water both as a CVP delivery and under their own right), name normalised as in valley3d_water_use_s2.py, exact or a
#      unique start of at least 12 characters; or
#  (b) it repeats, for the same owner, source and year, the exact amount of another right of that owner (one statement filed
#      twice); the first by application number is kept in the total.
import re  # noqa: E402
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
    for a_, b_ in ABBR:
        s = re.sub(a_, b_, s)
    return re.sub(r"\s+", " ", s).strip()


W2 = json.load(open(os.path.join(REPO, "figures", "valley_3d", "data", "water_use2.json"), encoding="utf-8"))
cvpn = {norm(n) for d in W2["districts"] for n in d.get("cvp_names", [])} | {norm(n) for n in W2["unmatched_cvp"]}
cvpn = {k for k in cvpn if len(k) >= 6}
def is_cvp(owner):
    k = norm(owner)
    return k in cvpn or any(len(c) >= 12 and k.startswith(c) for c in cvpn)
for r_ in rights:
    r_["cvp_owner"] = bool(is_cvp(r_["owner"]))
seen = {}
for r_ in sorted(rights, key=lambda x: x["appl"]):
    for y_, v_ in r_["years"].items():
        key_ = (norm(r_["owner"]), str(r_["source"]).lower(), y_, v_["af"])
        if v_["af"] > 0 and key_ in seen:
            v_["dup_of"] = seen[key_]
        elif v_["af"] > 0:
            seen[key_] = r_["appl"]
doc["rights"] = rights
json.dump(doc, open(OUT, "w", encoding="utf-8"), separators=(",", ":"))
for y_ in (2010, 2015, 2019, 2021, 2023):
    t_all = sum(r_["years"].get(y_, {}).get("af", 0) for r_ in rights)
    t_cvp = sum(r_["years"].get(y_, {}).get("af", 0) for r_ in rights if r_["cvp_owner"])
    t_dup = sum(r_["years"].get(y_, {}).get("af", 0) for r_ in rights if not r_["cvp_owner"] and r_["years"].get(y_, {}).get("dup_of"))
    print(f"{y_}: drawn {t_all * 1233.48e-9:.2f} km3 | CVP-contractor owners {t_cvp * 1233.48e-9:.2f} | repeated amounts {t_dup * 1233.48e-9:.2f} | in total {(t_all - t_cvp - t_dup) * 1233.48e-9:.2f}")
print("rights with a CVP-contractor owner:", sum(r_["cvp_owner"] for r_ in rights), "e.g.", sorted({r_["owner"] for r_ in rights if r_["cvp_owner"]})[:12])
