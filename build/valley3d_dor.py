"""Degree of regulation (DOR) for the 3D page: per rim river and per dam. Page data only; no test.

DOR (%) = reservoir storage upstream / mean yearly natural flow x 100: how many per cent of a year's natural flow the dams can hold.

Per river (the page's 14 rim gauges):
  storage upstream of the gauge, by date: GAGES-II HydroMod_Dams (enhanced NID 2009), pre1940_STOR ... pre1990_STOR and STOR_NID_2009,
  in megaliters per km2 of the gauge basin, times DRAIN_SQKM (GAGES-II basinid); 1 ML = 1,000 m3. STOR_NOR_2009 (normal storage) is kept beside it.
  natural flow: mean of the decade means of CDEC full natural flow already on the page (data/case/case.json, rivers[].decades[].natural, km3 a year).
  flow before and after: natural (that mean) and flow measured below the dam in the 2010s (case.json), km3 a year.
Per dam: Global Dam Watch v1.0 barriers DOR_PC (storage over the long-term mean flow at the dam, modelled by GDW) and DIS_AVG_LS,
  matched to the page's dams (data/meta.json) by name and year built.
Writes figures/valley_3d/data/case/dor.json.
"""
import json, os, re
import pandas as pd
import geopandas as gpd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PAGE = os.path.join(REPO, "figures", "valley_3d", "data")
G2 = "Z:/Data/california_central_valley/gagesii_basinchar/csv"
GDW = "D:/Data/grand/gdw_v1_0/GDW_v1_0_shp/GDW_barriers_v1_0.shp"

case = json.load(open(os.path.join(PAGE, "case", "case.json"), encoding="utf-8"))
meta = json.load(open(os.path.join(PAGE, "meta.json"), encoding="utf-8"))
dams = pd.read_csv(os.path.join(G2, "conterm_hydromod_dams.txt"), dtype={"STAID": str})
bas = pd.read_csv(os.path.join(G2, "conterm_basinid.txt"), dtype={"STAID": str}, encoding="latin-1")
STOR = [("before 1940", "pre1940_STOR"), ("before 1950", "pre1950_STOR"), ("before 1960", "pre1960_STOR"), ("before 1970", "pre1970_STOR"),
        ("before 1980", "pre1980_STOR"), ("before 1990", "pre1990_STOR"), ("2009", "STOR_NID_2009")]

rivers = []
for r in case["rivers"]:
    sid = str(r["site"]); d = dams[dams.STAID == sid]; b = bas[bas.STAID == sid]
    nat = [v["natural"] for v in r["decades"].values() if v.get("natural") is not None]
    if d.empty or b.empty or not nat:
        print("skip", r["river"], sid); continue
    A = float(b.DRAIN_SQKM.iloc[0]); Q = sum(nat) / len(nat)
    km3 = lambda col: float(d[col].iloc[0]) * A * 1e-6          # ML/km2 x km2 = ML; 1 ML = 1e-6 km3
    series = [{"when": w, "storage_km3": round(km3(c), 3), "dor_pct": round(100 * km3(c) / Q, 1)} for w, c in STOR]
    m10 = r["decades"].get("2010", {}).get("measured")
    rivers.append({"river": r["river"], "site": sid, "gauge": b.STANAME.iloc[0], "basin_km2": round(A, 1),
                   "natural_km3_yr": round(Q, 3), "natural_decades": len(nat), "measured_2010s_km3_yr": m10,
                   "dor_by_date": series, "storage_normal_2009_km3": round(km3("STOR_NOR_2009"), 3),
                   "dor_normal_2009_pct": round(100 * km3("STOR_NOR_2009") / Q, 1), "dams_2009": int(d.NDAMS_2009.iloc[0])})
    print(f"{r['river']:12s} {sid} A {A:9.1f} km2  natural {Q:6.3f} km3/yr  DOR " + "  ".join(f"{s['when']} {s['dor_pct']:.1f}" for s in series)
          + f"  | normal 2009 {rivers[-1]['dor_normal_2009_pct']:.1f}  | 2010s below dam {m10}")

g = gpd.read_file(GDW, bbox=(-125.0, 34.0, -117.0, 42.1))   # lon/lat box around the page grid
key = lambda s: re.sub(r"[^a-z0-9]", "", str(s).lower().replace(" dam", ""))
gk = {}
for _, x in g.iterrows(): gk.setdefault((key(x.DAM_NAME), int(x.YEAR_DAM)), x)
out_d, miss = [], 0
for d in meta["dams"]:
    x = gk.get((key(d["name"]), int(d["built"] or -1)))
    if x is None or x.DOR_PC is None or x.DOR_PC < 0: miss += 1; continue
    out_d.append({"name": d["name"], "built": d["built"], "dor_pct": round(float(x.DOR_PC), 1), "mean_flow_m3s": round(float(x.DIS_AVG_LS) / 1000, 1)})
print(f"dams with a GDW DOR: {len(out_d)} of {len(meta['dams'])}; not matched or no value: {miss}")
json.dump({"SOURCE": __doc__.strip(), "rivers": rivers, "dams": out_d}, open(os.path.join(PAGE, "case", "dor.json"), "w", encoding="utf-8"), ensure_ascii=False)
print("wrote", os.path.join(PAGE, "case", "dor.json"))
