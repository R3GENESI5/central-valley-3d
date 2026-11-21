"""Low-level winds for the 3D page's wind layer: how the wind runs across the valley on storm days and on the afternoon after,
winters 1959-60 to 1979-80 against 1999-2000 to 2021-22 (Ali, 2026-11-18: "the purpose of this was to show how the wind runs
across the valley" ... "Yes, download"). Page data, not a test.

Days, eras, grid and floor cells: those of src/landscape/valley3d_wind.py (imported, not copied): storm day = offshore inflow
>= 150 kg/m/s (corridor rule, era5_threshold), all four steps; afternoon after = the 00 UTC step of the first day without inflow
after a storm day; November to March; water years 1960-1980 and 2000-2022; the page's 5 km grid.
Winds: ERA5 0.25 degree u, v at 10 m (single levels) and at 925 and 850 hPa (pressure levels), from
Z:/Data/california_central_valley/era5_lowwind_em (src/level0/fetch_lowwind_em.py, Earthmover's public copy on AWS; equal to the
CDS files, raw/tests/2026-11-18_lowwind-earthmover-check; the files hold only the storm-day and afternoon-after steps). Where 925 hPa lies below the ground (the Sierra
above about 700 m) ERA5 extrapolates; those cells are kept on the grid and the page says so.
ERA5's land cover and roughness are one fixed map for all years, so a 10 m wind change here comes from the weather and the
analysed air, not from a change in roughness.
Writes figures/valley_3d/data/wind/wind_low.bin (float32; era, composite, level [10 m, 925, 850], component [u, v], row, column)
and wind_low_meta.json. Prints the floor-mean wind (speed, direction it blows from, along-valley component toward 330 degrees)
per era, composite and level.
Run  python src/landscape/valley3d_lowwind.py
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import xarray as xr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import valley3d_wind as W  # noqa: E402  (grid, floor, eras, day rule; its main() is not run)

Z = W.Z
LW = os.path.join(Z, "era5_lowwind_em")   # Earthmover AWS copy, equal to CDS (raw/tests/2026-11-18_lowwind-earthmover-check)
OUT = W.OUT
LEVELS = ["10m", "925", "850"]
AXIS = 330.0   # degrees: the valley's long axis, pointing north-north-west


def winter(kind, wy):
    parts = []
    for y, months in ((wy - 1, [10, 11, 12]), (wy, [1, 2, 3])):
        d = xr.open_dataset(os.path.join(LW, f"era5_{kind}_wind_{y}.nc"))
        t = "valid_time" if "valid_time" in d.dims else "time"
        d = d.rename({t: "time"}) if t != "time" else d
        parts.append(d.sel(time=d.time.dt.month.isin(months)).load())
    return xr.concat(parts, "time")


def days_of(wy):
    """Storm days and afternoons after, exactly as valley3d_wind.py (offshore inflow from era5_threshold)."""
    b = np.radians(W.BEARING)
    parts = []
    for y in (wy - 1, wy):
        ds = xr.open_dataset(os.path.join(Z, "era5_threshold", f"era5_ivt_t2_td2_{y}.nc"))[["viwve", "viwvn"]]
        parts.append(ds.sel(time=ds.time.dt.month.isin([10, 11, 12] if y == wy - 1 else [1, 2, 3, 4])).load())
    T = xr.concat(parts, "time")
    LO2, LA2 = np.meshgrid(T.longitude.values, T.latitude.values)
    al, cr = W.along_cross(LA2, LO2)
    kb = np.where((al >= 0) & (al < 3 * W.BIN_KM) & (np.abs(cr) <= W.HALF_KM), np.floor(al / W.BIN_KM).astype(int), -1)
    fa = T.viwve.values * np.sin(b) + T.viwvn.values * np.cos(b)
    f0 = np.mean([fa[:, kb == q].mean(1) for q in range(3)], axis=0)
    t = pd.DatetimeIndex(T.time.values)
    pdy = W.pday(t)
    g = pd.DataFrame({"f": f0, "d": pdy}).groupby("d").f.agg(["mean", "size"])
    g = g[g["size"] == 4]
    g = g[g.index.month.isin([11, 12, 1, 2, 3]) & (np.where(g.index.month >= 10, g.index.year + 1, g.index.year) == wy)]
    inflow = set(g.index[g["mean"] >= 150])
    after = {d + pd.Timedelta(days=1) for d in inflow if (d + pd.Timedelta(days=1)) in g.index and (d + pd.Timedelta(days=1)) not in inflow}
    storm_t = [x for x, d in zip(t, pdy) if d in inflow]
    after_t = [x for x, d in zip(t, pdy) if d in after and x.hour == 0]
    return {"storm": storm_t, "after": after_t}


def main():
    if os.path.exists(os.path.join(OUT, "wind_low.bin")):
        raise SystemExit(f"{OUT}/wind_low.bin exists; not overwritten")
    S = np.zeros((2, 2, 3, 2, W.WY, W.WX)); N = np.zeros((2, 2), int); miss = 0
    for ei, (era, wys) in enumerate(W.ERAS.items()):
        for wy in wys:
            D = days_of(wy)
            sl, pl = winter("sl", wy), winter("pl", wy)
            st, pt = pd.DatetimeIndex(sl.time.values), pd.DatetimeIndex(pl.time.values)
            la, lo = sl.latitude.values, sl.longitude.values
            lo = np.where(lo > 180, lo - 360, lo)
            for ci, c in enumerate(W.COMPS):
                ts = [x for x in D[c] if x in st and x in pt]
                miss += len(D[c]) - len(ts)
                if not ts:
                    continue
                i_s, i_p = [st.get_loc(x) for x in ts], [pt.get_loc(x) for x in ts]
                fields = [sl.u10.values[i_s].mean(0), sl.v10.values[i_s].mean(0)]
                for lev in (925, 850):
                    p = pl.sel(pressure_level=lev)
                    fields += [p.u.values[i_p].mean(0), p.v.values[i_p].mean(0)]
                for k, f in enumerate(fields):
                    S[ei, ci, k // 2, k % 2] += W.to_grid(f, la, lo) * len(ts)
                N[ei, ci] += len(ts)
            print(era, wy, "storm steps", len(D["storm"]), "afternoons", len(D["after"]), flush=True)
    A = (S / np.maximum(1, N)[:, :, None, None, None, None]).astype(np.float32)
    A.tofile(os.path.join(OUT, "wind_low.bin"))
    ax = np.radians(AXIS)
    res = {}
    for ei, era in enumerate(W.ERAS):
        for ci, c in enumerate(W.COMPS):
            for li, lev in enumerate(LEVELS):
                u, v = A[ei, ci, li, 0][W.FLOOR].mean(), A[ei, ci, li, 1][W.FLOOR].mean()
                spd = np.hypot(A[ei, ci, li, 0], A[ei, ci, li, 1])[W.FLOOR].mean()
                frm = (np.degrees(np.arctan2(-u, -v)) + 360) % 360
                along = u * np.sin(ax) + v * np.cos(ax)
                res[f"{era}_{c}_{lev}"] = [round(float(spd), 2), round(float(frm)), round(float(along), 2)]
                print(f"{era:6s} {c:5s} {lev:4s} floor mean speed {spd:5.2f} m/s | mean wind from {frm:5.0f} deg | along-valley (toward {AXIS:.0f} deg) {along:+5.2f} m/s")
    steps = {era: {c: int(N[ei, ci]) for ci, c in enumerate(W.COMPS)} for ei, era in enumerate(W.ERAS)}
    print("six-hourly steps used:", steps, "; steps without a wind file match:", miss)
    json.dump({"SOURCE": __doc__.strip(), "order": {"era": list(W.ERAS), "comp": W.COMPS, "level": LEVELS, "component": ["u", "v"]},
               "steps": steps, "floor": res, "floor_key": "speed m/s, direction the mean wind blows from (deg), along-valley component toward 330 deg (m/s)"},
              open(os.path.join(OUT, "wind_low_meta.json"), "w", encoding="utf-8"), ensure_ascii=False)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
