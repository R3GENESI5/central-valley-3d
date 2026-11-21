"""Stage 1 data for the "Where the water goes" support layer of the 3D Central Valley page (Ali, 2026-11-20: water
contracts layer, "plan first ... optional or additive layer, not directly part of the diagnostic"; plan agreed, "Yes,
stage 1"). Only data already on Z:, nothing downloaded.

Writes figures/valley_3d/data/water_use.json:
  rivers  full natural flow per rim river, water years (Oct-Sep) in acre-feet, CDEC FNF monthly, each drawn at the dam
          its gauge sits at or just below (anchor named in the file)
  canals  CVP deliveries per canal system, calendar years 1998-2025, 12 months + printed annual in acre-feet, from the
          USBR CVO monthly tables 22 (Friant-Kern), 23 (Madera and Millerton), 25 (Delta-Mendota), 26 (San Luis and
          Cross Valley), 27 (Tehama-Colusa), 28 (Sacramento River); the system total is the row labelled Total with the
          largest annual; a year counts only when its 12 months sum to the printed annual within 2 AF
  delta   Dayflow daily, from 1929-10-01: OUT (net Delta outflow), EXP (total exports), SWP (Banks) and CVP (Jones) exports
          (from 1956; earlier files give the total only), TOT (inflow), cfs
SOURCE  Z:/Data/california_central_valley/water_balance/cdec_fnf_monthly.csv; .../imports/usbr_cvo/raw (USBR CVO
        tables, fetched 2026-10-30); .../dayflow/dayflow-results-*.csv (DWR Dayflow).
WHEN    2026-11-20. STATUS derived (parsing, unit sums; nothing fitted). Row readers from src/level0/parse_cvp_deliveries.py
        (imported, not edited).
Run     python -P src/landscape/valley3d_water_use.py
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "src", "level0"))
import parse_cvp_deliveries as P  # noqa: E402  (rows_txt, rows_pdf)

Z = r"Z:\Data\california_central_valley"
OUT = os.path.join(REPO, "figures", "valley_3d", "data", "water_use.json")

# CDEC FNF station -> (river, page dam it is drawn at)
FNF = {"SBB": ("Sacramento River at Bend Bridge", "Shasta"), "FTO": ("Feather River at Oroville", "Oroville"),
       "YRS": ("Yuba River near Smartsville", "Harry L. Englebright Dam"), "AMF": ("American River below Folsom", "Folsom"),
       "CSN": ("Cosumnes River at Michigan Bar", None), "MKM": ("Mokelumne River inflow to Pardee", "Pardee"),
       "SNS": ("Stanislaus River below Goodwin", "New Melones"), "TLG": ("Tuolumne River below La Grange", "Don Pedro"),
       "MRC": ("Merced River near Merced Falls", "Exchequer Main"), "SBF": ("San Joaquin River below Friant", "Friant"),
       "KGF": ("Kings River below Pine Flat", "Pine Flat Dam"), "KWT": ("Kaweah River below Terminus", "Terminus Dam"),
       "SCC": ("Tule River below Success", "Success Dam"), "KRI": ("Kern River below Isabella", "Isabella Dam")}
CANALS = {22: ("fkc", "Friant-Kern Canal", "Friant-Kern Canal"), 23: ("madera", "Madera Canal and Millerton Lake", "Madera Canal"),
          25: ("dmc", "Delta-Mendota Canal", "Delta-Mendota Canal"),
          26: ("sanluis", "San Luis and Cross Valley Canals", "Governor Edmund G Brown California Aqueduct"),
          27: ("tcc", "Tehama-Colusa Canal", None), 28: ("sacriver", "Sacramento River (settlement and other contractors)", None)}

# rivers
f = pd.read_csv(os.path.join(Z, "water_balance", "cdec_fnf_monthly.csv"), parse_dates=["month"])
f["wy"] = f["month"].dt.year + (f["month"].dt.month >= 10)
rivers = {}
for st, (name, dam) in FNF.items():
    g = f[f.station == st].groupby("wy")["fnf_af"].agg(["sum", "count"])
    g = g[g["count"] == 12]
    rivers[st] = {"name": name, "dam": dam, "wy_af": {int(k): round(float(v)) for k, v in g["sum"].items()}}

# canals
def tfile(t, y):
    yy = y % 100
    if y <= 2000:
        return os.path.join(Z, "imports", "usbr_cvo", "raw", f"tab{t}{yy:02d}.txt")
    if y <= 2010:
        return os.path.join(Z, "imports", "usbr_cvo", "raw", f"tab{t}{yy:02d}.pdf")
    return os.path.join(Z, "imports", "usbr_cvo", "raw", f"table_{t}_{y}.pdf")


canals, log = {}, []
for t, (key, name, line) in CANALS.items():
    canals[key] = {"name": name, "line": line, "table": t, "years": {}}
    for y in range(1998, 2026):
        fn = tfile(t, y)
        if not os.path.exists(fn):
            log.append(f"{key} {y}: no file"); continue
        rows = P.rows_txt(fn) if fn.endswith(".txt") else P.rows_pdf(fn)
        tots = [r for r in rows if r[0].strip().lower().startswith("total") and r[2] is not None]
        if not tots:
            log.append(f"{key} {y}: no Total row"); continue
        lab, m, tot = max(tots, key=lambda r: r[2])
        m = [0.0 if v is None else v for v in m]
        if abs(sum(m) - tot) > 2:
            log.append(f"{key} {y}: months {sum(m):.0f} vs printed {tot:.0f}"); continue
        canals[key]["years"][y] = {"af": round(tot), "m": [round(v) for v in m], "row": lab.strip()}

# Delta (Dayflow)
fr = []
for fn in sorted(glob.glob(os.path.join(Z, "dayflow", "dayflow-results-*.csv"))):
    d = pd.read_csv(fn)
    d = d.rename(columns={"OUT1": "OUT", "TOT1": "TOT", "OUT2": "OUT", "TOT2": "TOT", "EXPOR": "EXP", "EXPORT": "EXP", "EXPORTS": "EXP"})
    d["date"] = pd.to_datetime(d["Date"])
    for c in ("SWP", "CVP", "EXP"):   # 1929-39 has no exports column (none then); 1940-55 gives total exports only
        if c not in d:
            d[c] = np.nan if c != "EXP" else 0.0
    fr.append(d[["date", "OUT", "SWP", "CVP", "EXP", "TOT"]])
dd = pd.concat(fr).drop_duplicates("date").sort_values("date").set_index("date")
full = pd.date_range(dd.index[0], dd.index[-1], freq="D")
dd = dd.reindex(full)
delta = {"start": str(full[0].date()), "unit": "cfs", "n": len(full),
         **{c.lower(): [None if pd.isna(v) else int(round(v)) for v in dd[c]] for c in ("OUT", "SWP", "CVP", "EXP", "TOT")}}

doc = {"SOURCE": "CDEC full natural flow (water_balance/cdec_fnf_monthly.csv); USBR CVO monthly delivery tables 22-28, 1998-2025 "
                 "(imports/usbr_cvo/raw); DWR Dayflow 1929-2025 (dayflow)", "WHEN": "2026-11-20",
       "STATUS": "derived (parsing, unit sums; nothing fitted)", "AF_M3": 1233.48, "CFS_M3S": 0.0283168,
       "rivers": rivers, "canals": canals, "delta": delta, "parse_log": log}
json.dump(doc, open(OUT, "w", encoding="utf-8"), separators=(",", ":"))

print("size", os.path.getsize(OUT) // 1024, "KB")
for st, r in rivers.items():
    ys = sorted(r["wy_af"]); print(f"FNF {st} {r['name']}: WY {ys[0]}-{ys[-1]} ({len(ys)}), WY2019 {r['wy_af'].get(2019, 'na'):,} AF, drawn at {r['dam']}")
for k, c in canals.items():
    ys = sorted(c["years"]); print(f"canal {k} ({c['name']}): {len(ys)} years {ys[0] if ys else '-'}-{ys[-1] if ys else '-'}, 2019 {c['years'].get(2019, {}).get('af', 'na')} AF, row '{c['years'].get(2019, {}).get('row', '')}'")
print("parse log:", len(log)); [print("  ", l) for l in log]
print("dayflow", delta["start"], "to", str(full[-1].date()), delta["n"], "days; missing OUT", sum(v is None for v in delta["out"]))
y19 = dd.loc["2019"]; k = 0.0283168 * 86400 / 1e9
print("2019 km3: OUT %.2f SWP %.2f CVP %.2f EXP %.2f TOT %.2f" % tuple(y19[c].sum() * k for c in ("OUT", "SWP", "CVP", "EXP", "TOT")))
