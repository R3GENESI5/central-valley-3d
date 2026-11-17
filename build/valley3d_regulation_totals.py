"""Regional totals of river regulation for the 3D page, by date: page data only, no test.

Region: the page grid (the Central Valley and the mountains around it).
Dates: before 1900, before 1940, before 1970, now (all dams with a build year).
  dams:            page dams (Global Dam Watch v1.0 in data/meta.json) inside the grid, built before the date.
  storage:         sum of their storage (GDW CAP_MCM as stored in meta.json), km3.
  streams dammed:  distinct river names of those dams (GDW RIVER field, tributary forks folded into their river).
  major rivers:    of the page's 14 gauged rim rivers, how many had any dam storage, storage of at least 10 % and at least 100 % of the yearly natural flow, upstream of the gauge by the date (GAGES-II pre1940_STOR, pre1970_STOR, STOR_NID_2009; no field before 1940).
  regional DOR:    storage upstream of the 14 gauges (GAGES-II) over their summed mean natural flow (CDEC, case.json), per cent.
Also: named rivers in the page's river network (data/rivers.json: NHD flowlines on the floor and lower foothills, Natural Earth beyond), and the
summed natural flow against the summed flow measured below the dams in the 2010s, for the rivers that have both.
Writes figures/valley_3d/data/case/totals.json.
"""
import json, os, re
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PAGE = os.path.join(REPO, "figures", "valley_3d", "data")
G2 = "Z:/Data/california_central_valley/gagesii_basinchar/csv"
meta = json.load(open(os.path.join(PAGE, "meta.json"), encoding="utf-8"))
case = json.load(open(os.path.join(PAGE, "case", "case.json"), encoding="utf-8"))
riv = json.load(open(os.path.join(PAGE, "rivers.json"), encoding="utf-8"))
TW, TH = meta["grid"]["tw"], meta["grid"]["th"]

def stream(n):
    n = (n or "").lower().strip()
    n = re.sub(r"^(north|south|middle|east|west|little) fork (of )?(the )?", "", n)
    n = re.sub(r"\b(north|south|middle|east|west) fork\b", "", n)
    return re.sub(r"\s+", " ", n).strip()

dams = [d for d in meta["dams"] if d.get("built") and 0 <= d["px"][0] < TW and 0 <= d["px"][1] < TH]
print(f"page dams inside the grid with a build year: {len(dams)} (of {len(meta['dams'])} in meta.json)")
dd = pd.read_csv(os.path.join(G2, "conterm_hydromod_dams.txt"), dtype={"STAID": str})
bas = pd.read_csv(os.path.join(G2, "conterm_basinid.txt"), dtype={"STAID": str}, encoding="latin-1")
R = []
for r in case["rivers"]:
    sid = str(r["site"]); d = dd[dd.STAID == sid]; b = bas[bas.STAID == sid]
    nat = [v["natural"] for v in r["decades"].values() if v.get("natural") is not None]
    A = float(b.DRAIN_SQKM.iloc[0])
    R.append({"river": r["river"], "Q": sum(nat) / len(nat), "m10": r["decades"].get("2010", {}).get("measured"),
              "s": {k: float(d[c].iloc[0]) * A * 1e-6 for k, c in [("1940", "pre1940_STOR"), ("1970", "pre1970_STOR"), ("now", "STOR_NID_2009")]}})
Qsum = sum(x["Q"] for x in R)
rows = []
for label, cut, g2 in [("before 1900", 1900, None), ("before 1940", 1940, "1940"), ("before 1970", 1970, "1970"), ("now", 9999, "now")]:
    ds = [d for d in dams if d["built"] < cut]
    row = {"when": label, "dams": len(ds), "storage_km3": round(sum(d.get("mcm") or 0 for d in ds) / 1000, 2),
           "streams_dammed": len({stream(d.get("river")) for d in ds if stream(d.get("river"))})}
    if g2:
        S = sum(x["s"][g2] for x in R)
        row.update({"major_rivers_dammed": sum(1 for x in R if x["s"][g2] > 0), "major_rivers_dor10": sum(1 for x in R if x["s"][g2] >= 0.1 * x["Q"]), "major_rivers_dor100": sum(1 for x in R if x["s"][g2] >= x["Q"]), "major_rivers": len(R), "storage_14_km3": round(S, 2), "regional_dor_pct": round(100 * S / Qsum, 1)})
    rows.append(row); print(row)
both = [x for x in R if x["m10"] is not None]
flow = {"rivers_with_both": len(both), "natural_km3_yr": round(sum(x["Q"] for x in both), 2), "below_dams_2010s_km3_yr": round(sum(x["m10"] for x in both), 2)}
flow["passed_pct"] = round(100 * flow["below_dams_2010s_km3_yr"] / flow["natural_km3_yr"], 1)
print("14 rim rivers, summed mean natural flow", round(Qsum, 2), "km3/yr;", flow)
names = {stream(l.get("n")) for l in riv["lines"] if l.get("n") and " river" in (l.get("n") or "").lower()}
names.discard("")
dammed_now = {stream(d.get("river")) for d in dams}
print(f"named rivers in the page network: {len(names)}; of these with a page dam on them (name match): {len(names & dammed_now)}")
json.dump({"SOURCE": __doc__.strip(), "by_date": rows, "flow_14": flow, "natural_14_km3_yr": round(Qsum, 2),
           "named_rivers_network": len(names), "named_rivers_with_dam": len(names & dammed_now)},
          open(os.path.join(PAGE, "case", "totals.json"), "w", encoding="utf-8"), ensure_ascii=False)
print("wrote", os.path.join(PAGE, "case", "totals.json"))
