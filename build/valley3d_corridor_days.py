"""The week after a winter rain, day by day, 1982-99 against 2000-22, at the corridor's stations, for the 3D page
(Ali, 2026-11-18: "your corridor strip units or change is so close in pre and post that it conveys nothing" ... "go").
Page data, not a test: the measured day-by-day values behind the slopes of raw/tests/2026-11-17_corridor-full-ladder.

Same stations, groups, days, rain rule and quantities as src/corridor/corridor_full_ladder.py, imported (station(), GRP, the
daily table, corridor_bc.dsr_of and buoy(); nothing edited). For each group and era (water years 1982-1999 and 2000-2022) and each
day since rain d = 1..8: the mean of the quantity over all days with that d. Shown as the rise from day 1 (value on day d minus
value on day 1, same era), so an instrument step between eras (airport dew-point sensors, 1990s) does not enter.
Quantities: afternoon temperature (13-15 PST mean), afternoon temperature minus dew point, fog hours a day (visibility below
1,000 m); the buoy's daily air temperature.
Range: the rise by day 7, era by era and its change, 2,000 resamples of water years (same seed as the ladder).
Writes figures/valley_3d/data/case/corridor_days.json; prints every value.
Run  python src/landscape/valley3d_corridor_days.py
"""
import json
import os
import sys

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "src", "corridor"))
import corridor_full_ladder as L  # noqa: E402
cb = L.cb

OUT = os.path.join(REPO, "figures", "valley_3d", "data", "case", "corridor_days.json")
ERAS = {"1982-99": np.arange(1982, 2000), "2000-22": np.arange(2000, 2023)}
DAYS = np.arange(1, 9)
Q = {"aT": "afternoon temperature (K)", "agap": "afternoon temperature minus dew point (K)", "fog": "fog hours a day"}


def curves(y, d, wy):
    """per water year and d: sums and counts; returns arrays [41 years, 8 days]."""
    k = np.isfinite(y) & (d >= 1) & (d <= 8) & (wy >= 1982) & (wy <= 2022)
    s = np.zeros((41, 8)); n = np.zeros((41, 8))
    np.add.at(s, ((wy[k] - 1982).astype(int), (d[k] - 1).astype(int)), y[k]); np.add.at(n, ((wy[k] - 1982).astype(int), (d[k] - 1).astype(int)), 1)
    return s, n


def era_curve(s, n, years):
    S, N = s[years - 1982].sum(0), n[years - 1982].sum(0)
    return np.where(N >= 20, S / np.maximum(N, 1), np.nan), N


def main():
    if os.path.exists(OUT):
        raise SystemExit(f"{OUT} exists; not overwritten")
    rd = lambda nm: pd.read_csv(f"{L.TAB}/{nm}.csv", index_col=0, parse_dates=True)
    P = rd("P"); cols = list(P.columns); days = P.index
    wy = np.where(days.month >= 10, days.year + 1, days.year)
    seg_p = {g: P[[cols[i] for i in b]].mean(1, skipna=False) for g, (_, b) in L.GRP.items()}; seg_p["buoy"] = P[cols[2]]
    dsr = {g: np.asarray(cb.dsr_of(p).values, float) for g, p in seg_p.items()}
    series = {}
    for g, (sids, _) in L.GRP.items():
        X = pd.concat([L.station(s) for s in sids]).groupby(level=0).mean().reindex(days)
        for q in Q:
            series[(g, q)] = (X[q].values.astype(float), dsr[g])
    try:
        Bu = cb.buoy("46026").reindex(days); series[("buoy", "aT")] = (Bu.air.values.astype(float), dsr["buoy"])
    except Exception as ex:
        print("buoy not read:", type(ex).__name__, ex)
    rng = np.random.default_rng(L.SEED)
    out = {"SOURCE": __doc__.strip(), "days": DAYS.tolist(), "eras": list(ERAS), "quantities": Q, "groups": {}}
    print("group | quantity | era | rise from day 1 on days 1..8 | days counted on d=1..8")
    for (g, q), (y, d) in series.items():
        s, n = curves(y, d, wy); rec = {}
        for e, yrs in ERAS.items():
            c, N = era_curve(s, n, yrs); rise = c - c[0]
            rec[e] = {"rise": [None if not np.isfinite(v) else round(float(v), 3) for v in rise], "n": N.astype(int).tolist()}
            print(f"{g} | {q} | {e} | " + " ".join("  nan" if not np.isfinite(v) else f"{v:+.2f}" for v in rise) + " | " + " ".join(str(int(x)) for x in N))
        def r7(yrs, idx):
            S, N = s[idx - 1982].sum(0), n[idx - 1982].sum(0)
            return (S[6] / N[6] - S[0] / N[0]) if N[6] >= 20 and N[0] >= 20 else np.nan
        e0, e1 = ERAS["1982-99"], ERAS["2000-22"]
        bs = np.array([r7(None, rng.choice(e1, len(e1))) - r7(None, rng.choice(e0, len(e0))) for _ in range(L.NB)]); bs = bs[np.isfinite(bs)]
        b0 = np.array([r7(None, rng.choice(e0, len(e0))) for _ in range(L.NB)]); b1 = np.array([r7(None, rng.choice(e1, len(e1))) for _ in range(L.NB)])
        v0, v1 = r7(None, e0), r7(None, e1)
        J = lambda x: None if x is None or not np.isfinite(x) else round(float(x), 3)   # JSON has no NaN: too few days is written as null
        rec["day7"] = {"1982-99": [J(v0), J(np.nanpercentile(b0, 2.5)), J(np.nanpercentile(b0, 97.5))],
                       "2000-22": [J(v1), J(np.nanpercentile(b1, 2.5)), J(np.nanpercentile(b1, 97.5))],
                       "change": [J(v1 - v0), J(np.percentile(bs, 2.5)) if len(bs) else None, J(np.percentile(bs, 97.5)) if len(bs) else None]}
        f2 = lambda x: "  n/a" if x is None else f"{x:+.2f}"; D = rec["day7"]
        print(f"{g} | {q} | rise by day 7: 1982-99 {f2(D['1982-99'][0])} [{f2(D['1982-99'][1])}, {f2(D['1982-99'][2])}] | 2000-22 {f2(D['2000-22'][0])} [{f2(D['2000-22'][1])}, {f2(D['2000-22'][2])}] | change {f2(D['change'][0])} [{f2(D['change'][1])}, {f2(D['change'][2])}]")
        out["groups"].setdefault(g, {})[q] = rec
    json.dump(out, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1, allow_nan=False)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
