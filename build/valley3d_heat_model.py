"""Ground heat on the valley's hottest days, modelled for today's land and for the 1850 land (Ali, 2026-11-16: "ground
heat, needs to be modelled, some way, some method"; "Yes, go" to the energy-balance design).

Method, fixed before the model was run against MODIS:
  Days      the valley-wide hottest days of every year 2003-2025 (lst/lst_meta.json, 10 a year, 230 days).
  Observed  MODIS Aqua daytime LST (about 13:30) averaged over those days per pixel (at least half the days good),
            reprojected to the page's 1 km mesh; valley-floor nodes only.
  Air       floor-mean air temperature, dew point and wind at 20-21 UTC on those days, 12 ISD valley stations.
  Sun       clear-sky shortwave at 13:30 local solar time on the mean day of year, latitude of each node:
            Rs = 0.75 * 1361 * cos(zenith); daily clear-sky mean from FAO-56 Ra * 0.75.
  Energy    (1 - a) Rs + e Ld - e s Ts^4 - G - LE - rho cp (Ts - Ta) / (c * ra) = 0, solved for Ts by bisection;
            Ld Brutsaert (1975); G a share of Rn by class; ra from the log profile at 10 m wind with Monin-Obukhov stability
            (Businger-Dyer), z0 and d from the class; one scalar c on ra fitted on today's land, the only fitted number.
            r1 (neutral ra) is kept in the session scratchpad: RMSE 10.0 C, dry classes +13 C too hot.
  LE        instantaneous LE = daily ET (mm/d) * 28.36 W/m2 per mm/d * (Rs at 13:30 / daily mean Rs), capped at
            0.95 (Rn - G).
  ET today  independent of LST: CVHM2 Farm Process crop coefficients (Jun-Sep mean of the codes in each 2019 group)
            times Hargreaves reference ET from PRISM monthly Tmax/Tmin 2003-2025 (floor mean, Jun-Sep); land with no
            summer water (idle, dry native, alkali, grass) takes the Jun-Sep ET of the US-Var dry grassland tower.
            Check run with SSEBop Jun-Sep ET (p4_now, derived from LST, so its fit is partly circular).
  ET 1850   the project's pre-1900 Jun-Sep evaporation layer (p4_then: survival-budget class rates, 1922-1950 forcing).
  Classes   albedo, roughness and ground-flux share by class are literature-typical constants listed in PARAM
            (assumptions, not measured here).
  Cross-check  remnant method: today's observed LST of wet-native, water and dry-native nodes painted onto the 1850
            classes (wet types take the wet-native value, open water the water value, the rest the dry-native value).
Writes figures/valley_3d/data/lst/model_{then,now}.png and obs_clim.png (C = 10 + value / 4, 0 = none, 1 km mesh)
and lst/model_meta.json.
WHEN 2026-11-16. STATUS derived (one fitted scalar; everything else fixed constants or held data).
Run  python src/landscape/valley3d_heat_model.py
"""
import base64
import datetime as dt
import glob
import json
import math
import os
import zipfile

import h5py
import numpy as np
import pandas as pd
import rasterio
from PIL import Image
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject, transform

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(REPO, "figures", "valley_3d", "data")
Z = r"Z:\Data\california_central_valley"
META = json.load(open(os.path.join(DATA, "meta.json"), encoding="utf-8"))
G = META["grid"]; MW, MH = G["mw"], G["mh"]
LSTM = json.load(open(os.path.join(DATA, "lst", "lst_meta.json"), encoding="utf-8"))
DAYS = sorted({d["date"] for y in LSTM["years"].values() for d in y["days"]})
print("hottest days", len(DAYS))

# ---- floor nodes, 2019 group and 1850 type at each 1 km node
fl = np.array(Image.open(os.path.join(DATA, "floor.png"))).astype(np.int32)
cell = fl[..., 0] + 256 * fl[..., 1] - 1; nat = fl[..., 2]
jj, ii = np.mgrid[0:MH, 0:MW]
tx = np.clip(np.round(ii * G["mesh_m"] / G["tex_m"]).astype(int), 0, G["tw"] - 1); ty = np.clip(np.round(jj * G["mesh_m"] / G["tex_m"]).astype(int), 0, G["th"] - 1)
ncell, ntype = cell[ty, tx], nat[ty, tx]
floor = ncell >= 0
b = base64.b64decode(META["land"]["2019"]); codes = np.array([(b[k >> 1] & 15) if k & 1 else (b[k >> 1] >> 4) for k in range(META["cells"])])
GROUPS = META["groups"]          # water, wet_native, dry_native, rice, orchard, field, pasture, idle, irr_unspec, urban
grp = np.where(floor, codes[np.clip(ncell, 0, None)], -1)
TYPES = META["types"]            # none, grassland, oak, riparian, marsh, water, alkali, floodplain, chaparral
X = G["x0"] + ii * G["mesh_m"]; Y = G["y1"] - jj * G["mesh_m"]
lon, lat = transform(G["crs"], "EPSG:4326", X[floor].tolist(), Y[floor].tolist()); lat = np.array(lat)
print("floor nodes", int(floor.sum()))

# ---- observed: MODIS LST mean over the hottest days
EPOCH = dt.datetime(2000, 1, 1); want = set(DAYS)
acc = n = None; src_T = None
for path in sorted(glob.glob(os.path.join(Z, "modis_lst", "ca_cv_myd11a1_*", "MYD11A1.061_1km_aid0001.nc"))):
    with h5py.File(path, "r") as f:
        la, lo = f["lat"][:], f["lon"][:]
        src_T = from_origin(lo[0] - (lo[1] - lo[0]) / 2, la[0] - (la[1] - la[0]) / 2, lo[1] - lo[0], la[0] - la[1])
        for k, sec in enumerate(f["time"][:]):
            d = (EPOCH + dt.timedelta(seconds=int(sec))).date().isoformat()
            if d not in want:
                continue
            lst = f["LST_Day_1km"][k].astype(np.float64); qc = f["QC_Day"][k]
            good = (lst > 0) & ((qc & 3) <= 1) & (((qc >> 6) & 3) <= 1)
            v = np.where(good, lst * 0.02 - 273.15, 0.0)
            acc = v if acc is None else acc + v; n = good.astype(int) if n is None else n + good
obs_src = np.where(n >= len(DAYS) / 2, acc / np.maximum(n, 1), np.nan).astype(np.float32)
dst_T = from_origin(G["x0"] - G["mesh_m"] / 2, G["y1"] + G["mesh_m"] / 2, G["mesh_m"], G["mesh_m"])
obs = np.full((MH, MW), np.nan, np.float32)
reproject(obs_src, obs, src_transform=src_T, src_crs="EPSG:4326", dst_transform=dst_T, dst_crs=G["crs"], resampling=Resampling.bilinear, src_nodata=np.nan, dst_nodata=np.nan)
print("observed floor mean LST", round(float(np.nanmean(obs[floor])), 2), "C")

# ---- air: ISD stations at 20-21 UTC on those days
def isd_val(s, scale=10.0):
    v = s.astype(str).str.split(",").str[0]
    v = pd.to_numeric(v.str.replace("+", "", regex=False), errors="coerce")
    return v.where(v.abs() < 9000) / scale
rows = []
for y in sorted({d[:4] for d in DAYS}):
    for f in glob.glob(os.path.join(Z, "isd_hourly", y, "*.csv")):
        d = pd.read_csv(f, usecols=["STATION", "DATE", "TMP", "DEW", "WND"], low_memory=False)
        t = pd.to_datetime(d.DATE)
        sel = t.dt.strftime("%Y-%m-%d").isin(want) & t.dt.hour.isin([20, 21])
        d = d[sel]
        if len(d):
            w = d.WND.astype(str).str.split(",")
            rows.append(pd.DataFrame({"st": d.STATION, "T": isd_val(d.TMP), "Td": isd_val(d.DEW), "U": pd.to_numeric(w.str[3], errors="coerce").where(lambda s: s < 900) / 10.0}))
air = pd.concat(rows)
Ta, Td, U = float(air["T"].mean()), float(air["Td"].mean()), float(max(1.0, air["U"].mean()))
print(f"air at 20-21 UTC: T {Ta:.2f} C, Td {Td:.2f} C, wind {U:.2f} m/s ({len(air)} reports, {air.st.nunique()} stations)")

# ---- sun
doy = float(np.mean([dt.date.fromisoformat(d).timetuple().tm_yday for d in DAYS]))
dec = 0.409 * math.sin(2 * math.pi * doy / 365 - 1.39); dr = 1 + 0.033 * math.cos(2 * math.pi * doy / 365)
phi = np.radians(lat); omega = math.radians(15 * 1.5)                     # 13:30 solar time
cosz = np.sin(phi) * math.sin(dec) + np.cos(phi) * math.cos(dec) * math.cos(omega)
Rs = 0.75 * 1361 * cosz
ws = np.arccos(-np.tan(phi) * math.tan(dec))
Ra = 24 * 60 / math.pi * 0.0820 * dr * (ws * np.sin(phi) * math.sin(dec) + np.cos(phi) * math.cos(dec) * np.sin(ws))   # MJ/m2/d
Rs_day = 0.75 * Ra * 1e6 / 86400                                                                                            # W/m2
k_inst = Rs / Rs_day
print(f"mean day {doy:.0f}; Rs 13:30 {Rs.mean():.0f} W/m2, daily mean {Rs_day.mean():.0f}, ratio {k_inst.mean():.2f}")

# ---- ET today (independent of LST): Kc x Hargreaves ETo
def prism_floor_mean(var, ym):
    p = os.path.join(Z, "prism_monthly", "raw", var, f"prism_{var}_us_25m_{ym}.zip")
    tif = [x for x in zipfile.ZipFile(p).namelist() if x.endswith((".tif", ".bil"))][0]
    with rasterio.open("/vsizip/" + p.replace("\\", "/") + "/" + tif) as s:
        a = np.full((MH, MW), np.nan, np.float32)
        reproject(rasterio.band(s, 1), a, dst_transform=dst_T, dst_crs=G["crs"], resampling=Resampling.average, src_nodata=s.nodata, dst_nodata=np.nan)
    return float(np.nanmean(a[floor]))
eto_m = []
for m in (6, 7, 8, 9):
    tmx = np.mean([prism_floor_mean("tmax", f"{y}{m:02d}") for y in range(2003, 2026)])
    tmn = np.mean([prism_floor_mean("tmin", f"{y}{m:02d}") for y in range(2003, 2026)])
    dmid = dt.date(2015, m, 15).timetuple().tm_yday
    d_ = 0.409 * math.sin(2 * math.pi * dmid / 365 - 1.39); r_ = 1 + 0.033 * math.cos(2 * math.pi * dmid / 365); p_ = math.radians(37.0)
    w_ = math.acos(-math.tan(p_) * math.tan(d_)); ra = 24 * 60 / math.pi * 0.0820 * r_ * (w_ * math.sin(p_) * math.sin(d_) + math.cos(p_) * math.cos(d_) * math.sin(w_))
    eto_m.append(0.0023 * (ra / 2.45) * ((tmx + tmn) / 2 + 17.8) * math.sqrt(max(0, tmx - tmn)))
ETo = float(np.mean(eto_m))
kc = pd.read_csv(os.path.join(Z, "cvhm2_landuse", "CVHM2_CropCoefficients.txt"), sep="\t").set_index("Crop")[["June", "July", "August", "September"]].mean(1)
GROUP_CODES = {0: [1], 1: [23], 2: [3], 3: [17], 4: [4, 10, 12, 16], 5: [6, 7, 9, 11, 14, 18, 21, 22, 24], 6: [5, 13, 20], 7: [8], 8: [19], 9: [2, 15]}
tw = pd.read_csv(os.path.join(Z, "ameriflux", "tower_monthly_et.csv"))
dry_et = float(tw[(tw.site == "US-Var") & tw.month.isin([6, 7, 8, 9])].et_mm.mean() / 30.4)       # mm/d
kc_g = {g: float(np.mean([kc[c] for c in cs])) for g, cs in GROUP_CODES.items()}
et_g = {g: kc_g[g] * ETo for g in GROUP_CODES}
for g in (2, 7):                                                                                     # dry native, idle: no summer water
    et_g[g] = dry_et
et_g[9] = 0.4 * ETo                                                                                  # urban: assumed 0.4 (mixed paved and watered)
print(f"Hargreaves ETo Jun-Sep floor mean {ETo:.2f} mm/d; US-Var dry grassland Jun-Sep {dry_et:.2f} mm/d")
print("ET today by group (mm/d):", {GROUPS[g]: round(v, 2) for g, v in et_g.items()})

def season_daily(path):
    with rasterio.open(path) as s:
        a = np.full((MH, MW), np.nan, np.float32)
        reproject(rasterio.band(s, 1), a, dst_transform=dst_T, dst_crs=G["crs"], resampling=Resampling.nearest, src_nodata=s.nodata, dst_nodata=np.nan)
    return np.where(a >= 0, a / 122.0, np.nan)
BA = r"C:\Users\PC\Projects\map-plates\data\central_valley\ba"
et_then = season_daily(os.path.join(BA, "p4_then_jun_sep_mm_1km.tif"))
et_ssebop = season_daily(os.path.join(BA, "p4_now_jun_sep_mm_1km.tif"))

# ---- class parameters (assumptions): albedo, z0 (m), ground-flux share of Rn
PARAM = {"water": (0.08, 0.001, 0.40), "marsh": (0.15, 0.15, 0.10), "riparian": (0.14, 1.0, 0.08), "oak": (0.17, 0.5, 0.10),
         "grass": (0.23, 0.03, 0.20), "alkali": (0.28, 0.02, 0.25), "floodplain": (0.20, 0.05, 0.15), "chaparral": (0.16, 0.4, 0.12),
         "rice": (0.15, 0.10, 0.10), "orchard": (0.17, 0.60, 0.08), "field": (0.20, 0.08, 0.12), "pasture": (0.21, 0.04, 0.12),
         "idle": (0.25, 0.02, 0.25), "urban": (0.16, 0.80, 0.20)}
NOW_CLS = {0: "water", 1: "marsh", 2: "grass", 3: "rice", 4: "orchard", 5: "field", 6: "pasture", 7: "idle", 8: "field", 9: "urban"}
THEN_CLS = {0: "grass", 1: "grass", 2: "oak", 3: "riparian", 4: "marsh", 5: "water", 6: "alkali", 7: "floodplain", 8: "chaparral"}

SIG, EPS, RHOCP, K = 5.67e-8, 0.97, 1.16e3, 0.41
ea = 0.6108 * math.exp(17.27 * Td / (Td + 237.3))                        # kPa
Ld = 1.24 * (ea * 10 / (Ta + 273.15)) ** (1 / 7) * SIG * (Ta + 273.15) ** 4

def ra_stable(z0, ts, c):
    """Aerodynamic resistance with Monin-Obukhov stability (Businger-Dyer, unstable and stable), z0h = 0.1 z0m."""
    h = z0 / 0.123; d = 0.67 * h; zm = max(10.0, h + 2.0); z = zm - d
    L, ra = -1e6, None
    for _ in range(30):
        x = max(1e-6, 1 - 16 * z / L) ** 0.25 if L < 0 else 1.0
        if L < 0:
            psm = 2 * math.log((1 + x) / 2) + math.log((1 + x * x) / 2) - 2 * math.atan(x) + math.pi / 2
            psh = 2 * math.log((1 + x * x) / 2)
        else:
            psm = psh = -5 * min(z / L, 1.0)
        us = K * U / max(0.05, math.log(z / z0) - psm)
        ra = (math.log(z / (0.1 * z0)) - psh) / (K * us)
        hflux = RHOCP * (ts - Ta) / (c * ra)
        Lnew = -RHOCP * (Ta + 273.15) * us ** 3 / (K * 9.81 * hflux) if abs(hflux) > 1e-3 else 1e6
        if abs(Lnew - L) < 0.01 * abs(L): L = Lnew; break
        L = Lnew
    return ra

def solve(alb, z0, gs, et_day, kinst, rs, c):
    lo, hi = Ta - 15, Ta + 60
    for _ in range(40):
        ts = (lo + hi) / 2
        rn = (1 - alb) * rs + EPS * Ld - EPS * SIG * (ts + 273.15) ** 4
        g = gs * rn; le = min(et_day * 28.36 * kinst, 0.95 * max(rn - g, 0))
        hflux = RHOCP * (ts - Ta) / (c * ra_stable(z0, ts, c))
        if rn - g - le - hflux > 0: lo = ts
        else: hi = ts
    return (lo + hi) / 2

idx = np.argwhere(floor)
def run(cls_of, et_of, c):
    out = np.full((MH, MW), np.nan, np.float32)
    for q, (j, i) in enumerate(idx):
        cl = cls_of(j, i); e = et_of(j, i)
        if cl is None or not np.isfinite(e):
            continue
        a, z0, gs = PARAM[cl]
        out[j, i] = solve(a, z0, gs, e, k_inst[q], Rs[q], c)
    return out

now_cls = lambda j, i: NOW_CLS.get(int(grp[j, i]))
now_et_B = lambda j, i: et_g.get(int(grp[j, i]), np.nan)
def skill(model):
    m = np.isfinite(model) & np.isfinite(obs) & floor
    dlt = model[m] - obs[m]
    return {"bias_c": round(float(dlt.mean()), 2), "rmse_c": round(float(np.sqrt((dlt ** 2).mean())), 2), "r": round(float(np.corrcoef(model[m], obs[m])[0, 1]), 3), "n": int(m.sum())}
# fit c on today's land with the independent ET (minimise floor-mean bias, then report spread)
lo_c, hi_c = 0.05, 3.0
for _ in range(14):
    c = (lo_c + hi_c) / 2
    bias = skill(run(now_cls, now_et_B, c))["bias_c"]
    if bias > 0: hi_c = c                        # model too hot -> more exchange -> smaller resistance
    else: lo_c = c
C = (lo_c + hi_c) / 2
now_B = run(now_cls, now_et_B, C)
now_A = run(now_cls, lambda j, i: et_ssebop[j, i], C)
then = run(lambda j, i: THEN_CLS.get(int(ntype[j, i])) if ntype[j, i] > 0 else None, lambda j, i: et_then[j, i], C)
by_grp = {GROUPS[g]: [round(float(np.nanmean(obs[(grp == g) & floor])), 1), round(float(np.nanmean(now_B[(grp == g) & floor])), 1)] for g in range(len(GROUPS)) if ((grp == g) & floor).sum() >= 20}
print(f"fitted c {C:.3f}")
print("skill today, independent ET:", skill(now_B), "| SSEBop ET (partly circular):", skill(now_A))
print("observed vs model by 2019 group (C):", by_grp)

# remnant cross-check
wet_v = float(np.nanmean(obs[(grp == 1) & floor])); dry_v = float(np.nanmean(obs[(grp == 2) & floor])); wat_v = float(np.nanmean(obs[(grp == 0) & floor]))
remn = np.where(ntype == 5, wat_v, np.where(np.isin(ntype, [3, 4, 7]), wet_v, np.where(ntype > 0, dry_v, np.nan)))
remn = np.where(floor, remn, np.nan)
both = floor & np.isfinite(then) & np.isfinite(obs)
res = {"floor_mean_c": {"observed_today": round(float(np.nanmean(obs[both])), 2), "model_today_independent_et": round(float(np.nanmean(now_B[both])), 2),
                        "model_1850": round(float(np.nanmean(then[both])), 2), "remnant_1850": round(float(np.nanmean(remn[both])), 2)},
       "change_model_1850_minus_model_today_c": round(float(np.nanmean(then[both] - now_B[both])), 2),
       "change_remnant_minus_observed_c": round(float(np.nanmean(remn[both] - obs[both])), 2),
       "remnant_values_c": {"wet_native": round(wet_v, 2), "dry_native": round(dry_v, 2), "water": round(wat_v, 2)}}
print("results:", json.dumps(res))

def enc(a):
    return np.where(np.isfinite(a), np.clip(np.round((a - 10) * 4), 1, 255), 0).astype(np.uint8)
OUT = os.path.join(DATA, "lst")
Image.fromarray(enc(np.where(floor, then, np.nan)), "L").save(os.path.join(OUT, "model_then.png"))
Image.fromarray(enc(np.where(floor, now_B, np.nan)), "L").save(os.path.join(OUT, "model_now.png"))
Image.fromarray(enc(np.where(floor, obs, np.nan)), "L").save(os.path.join(OUT, "obs_clim.png"))
json.dump({"SOURCE": __doc__.split("WHEN")[0].strip(), "WHEN": "2026-11-16", "STATUS": "derived", "encoding": "C = 10 + value / 4; 0 = none",
           "days": len(DAYS), "air": {"T_c": round(Ta, 2), "Td_c": round(Td, 2), "wind_ms": round(U, 2), "reports": int(len(air))},
           "sun": {"doy": round(doy), "rs_1330_wm2": round(float(Rs.mean())), "rs_daily_wm2": round(float(Rs_day.mean()))},
           "eto_mm_d": round(ETo, 2), "dry_tower_et_mm_d": round(dry_et, 2), "et_today_by_group_mm_d": {GROUPS[g]: round(v, 2) for g, v in et_g.items()},
           "param": PARAM, "fitted_c": round(C, 3), "skill_independent_et": skill(now_B), "skill_ssebop_et_partly_circular": skill(now_A),
           "observed_vs_model_by_group_c": by_grp, **res}, open(os.path.join(OUT, "model_meta.json"), "w", encoding="utf-8"), indent=1)
