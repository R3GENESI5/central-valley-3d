"""2 m air temperature on storm days and on the afternoon after, winters 1959-60 to 1979-80 against 1999-2000 to 2021-22, for the
3D page's wind layer (Ali, 2026-11-18: "the wind being shown can carry the moisture ... and temperature can come as a shade of red").
Page data, not a test. Days, eras, grid and floor cells: those of src/landscape/valley3d_wind.py, through days_of() in
src/landscape/valley3d_lowwind.py (both imported, not copied, not edited). Field: ERA5 t2m (0.25 degree, era5_threshold), deg C.
Writes figures/valley_3d/data/wind/wind_t2m.bin (float32; era, composite, row, column) and wind_t2m_meta.json; prints the
floor means with water-year resample ranges.
Run  python src/landscape/valley3d_wind_t2m.py
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import xarray as xr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import valley3d_wind as W  # noqa: E402
from valley3d_lowwind import days_of  # noqa: E402

OUT = W.OUT


def main():
    if os.path.exists(os.path.join(OUT, "wind_t2m.bin")):
        raise SystemExit(f"{OUT}/wind_t2m.bin exists; not overwritten")
    S = np.zeros((2, 2, W.WY, W.WX)); N = np.zeros((2, 2), int); YR = []
    for ei, (era, wys) in enumerate(W.ERAS.items()):
        for wy in wys:
            D = days_of(wy)
            parts = []
            for y in (wy - 1, wy):
                ds = xr.open_dataset(os.path.join(W.Z, "era5_threshold", f"era5_ivt_t2_td2_{y}.nc"))[["t2m"]]
                parts.append(ds.sel(time=ds.time.dt.month.isin([10, 11, 12] if y == wy - 1 else [1, 2, 3, 4])).load())
            T = xr.concat(parts, "time"); t = pd.DatetimeIndex(T.time.values)
            for ci, c in enumerate(W.COMPS):
                idx = [t.get_loc(x) for x in D[c] if x in t]
                if not idx:
                    YR.append((era, wy, c, 0, 0.0)); continue
                g = W.to_grid(T.t2m.values[idx].mean(0) - 273.15, T.latitude.values, T.longitude.values)
                S[ei, ci] += g * len(idx); N[ei, ci] += len(idx)
                YR.append((era, wy, c, len(idx), float(np.nanmean(g[W.FLOOR])) * len(idx)))
            print(era, wy, flush=True)
    A = (S / np.maximum(1, N)[:, :, None, None]).astype(np.float32)
    A.tofile(os.path.join(OUT, "wind_t2m.bin"))
    Y = pd.DataFrame(YR, columns=["era", "wy", "comp", "n", "t2m"])
    rng = np.random.default_rng(W.SEED); res = {}
    for c in W.COMPS:
        y0, y1 = Y[(Y.era == "before") & (Y.comp == c)], Y[(Y.era == "now") & (Y.comp == c)]
        r = lambda d: d.t2m.sum() / d.n.sum()
        bs = [r(y1.iloc[rng.integers(0, len(y1), len(y1))]) - r(y0.iloc[rng.integers(0, len(y0), len(y0))]) for _ in range(W.NBOOT)]
        res[c] = [round(r(y0), 2), round(r(y1), 2), round(r(y1) - r(y0), 2), round(float(np.quantile(bs, .025)), 2), round(float(np.quantile(bs, .975)), 2)]
        print(f"{c:5s} floor 2 m temperature (deg C) before {res[c][0]} | now {res[c][1]} | change {res[c][2]:+} [{res[c][3]:+}, {res[c][4]:+}]")
    steps = {era: {c: int(N[ei, ci]) for ci, c in enumerate(W.COMPS)} for ei, era in enumerate(W.ERAS)}
    print("six-hourly steps used:", steps)
    json.dump({"SOURCE": __doc__.strip(), "order": {"era": list(W.ERAS), "comp": W.COMPS}, "steps": steps, "floor": res},
              open(os.path.join(OUT, "wind_t2m_meta.json"), "w", encoding="utf-8"), ensure_ascii=False)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
