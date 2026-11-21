"""Ground heat on the 1850 land, r3: groundwater-fed evaporation and three seasons (Ali, 2026-11-17: "Yes, rebuild",
after the r2 result traced to 1850 grass and alkali scrub given summer evaporation with no water table under them).

Rule, fixed before any r3 number was produced:
  Seasons   summer Jun-Sep, spring Mar-May, winter Nov-Feb (Nov-Dec counted with the following January).
  Days      per season-year 2003-2025, the 10 valley-wide hottest days by the r2 rule (valley3d_lst_peaks.py): MODIS
            Aqua MYD11A1 day LST, good pixel = QC bits 0-1 <= 1 and error bits 6-7 <= 1, candidate day >= 70 % of floor
            pixels good, rank by floor-mean LST. Summer reuses lst_meta.json (same rule, same days as r2).
  Energy    r2 unchanged (valley3d_heat_model.py): (1 - a) Rs + e Ld - e s Ts^4 - G - LE - H = 0, H with Monin-Obukhov
            stability, one scalar c on ra fitted per season on today's land, then observed = a + b * model fitted per
            season on today's land and applied to the 1850 run. Vectorised here; summer with groundwater off must
            reproduce r2 (c 2.341, RMSE 8.78 C) before r3 numbers are read.
  ET today  CVHM2 Kc (season months) x Hargreaves ETo (PRISM 2003-2025, season months); dry native and idle take the
            US-Var tower ET of the season months; urban 0.4 ETo.
  ET 1850   plate ba_04 layer p4_then_<season> (survival-budget class rates), then for classes that draw on a water table
            ET = max(p4 rate, Kc_native x ETo x max(0, 1 - d / d_ext)), the linear extinction-depth form of MODFLOW EVT
            (the form CVHM uses). Kc_native = CVHM2 native vegetation code 3 (1.0).
            Classes and extinction depths (assumptions): grassland, alkali scrub, floodplain 3 m (sensitivity 2 and 4 m);
            valley oak 6 m (sensitivity 4 and 8 m). Marsh, riparian, water keep the p4 rates.
  Water table 1850  depth d from the USGS Water-Supply Paper wells of 1905-14 on the floor (p2_wsp_wells.gpkg):
            measured depth, flowing wells set to 1 m; inverse-distance (power 2) mean of the 8 nearest wells within
            10 km of a node, no groundwater term beyond 10 km. 1905-14 is after the first pumping and canals, so these
            depths are the deep side of 1850 and the groundwater term is a lower bound.
  Remnant   today's LST of native patches painted on the 1850 classes, reported as a drained-state number: today's
            native ground sits over a floor water table at a median depth of 22 m (p2_stats.json).
Writes figures/valley_3d/data/lst/r3_<season>_{then,now,obs}.png (C = 10 + value / 4, 0 = none) and adds the season
to lst/model_r3_meta.json. r2 files are left as they are.
WHEN 2026-11-17. STATUS derived (two fitted numbers per season, a and b, plus c; constants listed as assumptions).
Run  python src/landscape/valley3d_heat_model_r3.py <summer|spring|winter> [--no-gw]
"""
import datetime as dt
import glob
import json
import math
import os
import sys
import zipfile

import fiona
import h5py
import numpy as np
import pandas as pd
import rasterio
from PIL import Image
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject, transform
from scipy.spatial import cKDTree

SEASON = sys.argv[1]; GW = "--no-gw" not in sys.argv
SEAS = {"summer": ([6, 7, 8, 9], "jun_sep", 122), "spring": ([3, 4, 5], "mar_may", 92), "winter": ([11, 12, 1, 2], "nov_feb", 120)}
MONTHS, P4, NDAY = SEAS[SEASON]
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(REPO, "figures", "valley_3d", "data"); OUT = os.path.join(DATA, "lst")
Z = r"Z:\Data\california_central_valley"
BA = r"C:\Users\PC\Projects\map-plates\data\central_valley\ba"
META = json.load(open(os.path.join(DATA, "meta.json"), encoding="utf-8"))
G = META["grid"]; MW, MH = G["mw"], G["mh"]
EPOCH = dt.datetime(2000, 1, 1)
CUBES = sorted(glob.glob(os.path.join(Z, "modis_lst", "ca_cv_myd11a1_*", "MYD11A1.061_1km_aid0001.nc")))

# ---- floor nodes
fl = np.array(Image.open(os.path.join(DATA, "floor.png"))).astype(np.int32)
cell = fl[..., 0] + 256 * fl[..., 1] - 1; nat = fl[..., 2]
jj, ii = np.mgrid[0:MH, 0:MW]
tx = np.clip(np.round(ii * G["mesh_m"] / G["tex_m"]).astype(int), 0, G["tw"] - 1); ty = np.clip(np.round(jj * G["mesh_m"] / G["tex_m"]).astype(int), 0, G["th"] - 1)
ncell, ntype = cell[ty, tx], nat[ty, tx]
floor = ncell >= 0
b = __import__("base64").b64decode(META["land"]["2019"]); codes = np.array([(b[k >> 1] & 15) if k & 1 else (b[k >> 1] >> 4) for k in range(META["cells"])])
GROUPS = META["groups"]; TYPES = META["types"]
grp = np.where(floor, codes[np.clip(ncell, 0, None)], -1)
X = G["x0"] + ii * G["mesh_m"]; Y = G["y1"] - jj * G["mesh_m"]
idx = np.argwhere(floor); FJ, FI = idx[:, 0], idx[:, 1]
lat = np.array(transform(G["crs"], "EPSG:4326", X[floor].tolist(), Y[floor].tolist())[1])
dst_T = from_origin(G["x0"] - G["mesh_m"] / 2, G["y1"] + G["mesh_m"] / 2, G["mesh_m"], G["mesh_m"])

# ---- days
def season_year(d):
    return d.year + 1 if d.month in (11, 12) else d.year
if SEASON == "summer":
    LSTM = json.load(open(os.path.join(OUT, "lst_meta.json"), encoding="utf-8"))
    DAYS = sorted({d["date"] for y in LSTM["years"].values() for d in y["days"]})
else:
    fm = None; cand = {}
    for path in CUBES:
        with h5py.File(path, "r") as f:
            if fm is None:
                la, lo = f["lat"][:], f["lon"][:]
                LO, LA = np.meshgrid(lo, la)
                xs, ys = transform("EPSG:4326", G["crs"], LO.ravel(), LA.ravel())
                px = ((np.array(xs) - G["x0"]) / G["tex_m"]).astype(int); py_ = ((G["y1"] - np.array(ys)) / G["tex_m"]).astype(int)
                ok = (px >= 0) & (px < G["tw"]) & (py_ >= 0) & (py_ < G["th"])
                m = np.zeros(px.size, bool); m[ok] = (fl[..., 0] + 256 * fl[..., 1])[py_[ok], px[ok]] > 0
                fm = m.reshape(LO.shape)
            for k, sec in enumerate(f["time"][:]):
                d = EPOCH + dt.timedelta(seconds=int(sec))
                if d.month not in MONTHS or not (2003 <= season_year(d) <= 2025):
                    continue
                lst = f["LST_Day_1km"][k]; qc = f["QC_Day"][k]
                gf = (lst > 0) & ((qc & 3) <= 1) & (((qc >> 6) & 3) <= 1) & fm
                if gf.sum() / fm.sum() >= 0.70:
                    cand.setdefault(season_year(d), []).append((d.date().isoformat(), float(lst[gf].mean() * 0.02 - 273.15)))
    DAYS = sorted(d for y in cand for d, _ in sorted(cand[y], key=lambda r: -r[1])[:10])
    print("season-years", len(cand), "candidate days per year (min, median)", min(map(len, cand.values())), int(np.median([len(v) for v in cand.values()])))
want = set(DAYS)
print(SEASON, "hottest days", len(DAYS), "| groundwater term", GW)

# ---- observed LST on those days
acc = n = None
for path in CUBES:
    with h5py.File(path, "r") as f:
        la, lo = f["lat"][:], f["lon"][:]
        src_T = from_origin(lo[0] - (lo[1] - lo[0]) / 2, la[0] - (la[1] - la[0]) / 2, lo[1] - lo[0], la[0] - la[1])
        for k, sec in enumerate(f["time"][:]):
            if (EPOCH + dt.timedelta(seconds=int(sec))).date().isoformat() not in want:
                continue
            lst = f["LST_Day_1km"][k].astype(np.float64); qc = f["QC_Day"][k]
            good = (lst > 0) & ((qc & 3) <= 1) & (((qc >> 6) & 3) <= 1)
            v = np.where(good, lst * 0.02 - 273.15, 0.0)
            acc = v if acc is None else acc + v; n = good.astype(int) if n is None else n + good
obs_src = np.where(n >= len(DAYS) / 2, acc / np.maximum(n, 1), np.nan).astype(np.float32)
obs = np.full((MH, MW), np.nan, np.float32)
reproject(obs_src, obs, src_transform=src_T, src_crs="EPSG:4326", dst_transform=dst_T, dst_crs=G["crs"], resampling=Resampling.bilinear, src_nodata=np.nan, dst_nodata=np.nan)
OBS = obs[FJ, FI]
print("observed floor mean LST", round(float(np.nanmean(OBS)), 2), "C")

# ---- air
def isd_val(s, scale=10.0):
    v = pd.to_numeric(s.astype(str).str.split(",").str[0].str.replace("+", "", regex=False), errors="coerce")
    return v.where(v.abs() < 9000) / scale
rows = []
for y in sorted({d[:4] for d in DAYS}):
    for f in glob.glob(os.path.join(Z, "isd_hourly", y, "*.csv")):
        d = pd.read_csv(f, usecols=["STATION", "DATE", "TMP", "DEW", "WND"], low_memory=False)
        t = pd.to_datetime(d.DATE); d = d[t.dt.strftime("%Y-%m-%d").isin(want) & t.dt.hour.isin([20, 21])]
        if len(d):
            w = d.WND.astype(str).str.split(",")
            rows.append(pd.DataFrame({"st": d.STATION, "T": isd_val(d.TMP), "Td": isd_val(d.DEW), "U": pd.to_numeric(w.str[3], errors="coerce").where(lambda s: s < 900) / 10.0}))
air = pd.concat(rows)
Ta, Td, U = float(air["T"].mean()), float(air["Td"].mean()), float(max(1.0, air["U"].mean()))
print(f"air at 20-21 UTC: T {Ta:.2f} C, Td {Td:.2f} C, wind {U:.2f} m/s ({len(air)} reports, {air.st.nunique()} stations)")

# ---- sun
doys = [dt.date.fromisoformat(d).timetuple().tm_yday for d in DAYS]
ang = np.angle(np.mean(np.exp(2j * np.pi * np.array(doys) / 365))); doy = (ang / (2 * np.pi) * 365) % 365   # circular mean (winter spans the new year)
dec = 0.409 * math.sin(2 * math.pi * doy / 365 - 1.39); dr = 1 + 0.033 * math.cos(2 * math.pi * doy / 365)
phi = np.radians(lat); omega = math.radians(15 * 1.5)
cosz = np.sin(phi) * math.sin(dec) + np.cos(phi) * math.cos(dec) * math.cos(omega)
Rs = 0.75 * 1361 * cosz
ws = np.arccos(-np.tan(phi) * math.tan(dec))
Ra = 24 * 60 / math.pi * 0.0820 * dr * (ws * np.sin(phi) * math.sin(dec) + np.cos(phi) * math.cos(dec) * np.sin(ws))
Rs_day = 0.75 * Ra * 1e6 / 86400; KINST = Rs / Rs_day
print(f"mean day {doy:.0f}; Rs 13:30 {Rs.mean():.0f} W/m2, daily mean {Rs_day.mean():.0f}, ratio {KINST.mean():.2f}")

# ---- ET today
def prism_floor_mean(var, ym):
    p = os.path.join(Z, "prism_monthly", "raw", var, f"prism_{var}_us_25m_{ym}.zip")
    tif = [x for x in zipfile.ZipFile(p).namelist() if x.endswith((".tif", ".bil"))][0]
    with rasterio.open("/vsizip/" + p.replace("\\", "/") + "/" + tif) as s:
        a = np.full((MH, MW), np.nan, np.float32)
        reproject(rasterio.band(s, 1), a, dst_transform=dst_T, dst_crs=G["crs"], resampling=Resampling.average, src_nodata=s.nodata, dst_nodata=np.nan)
    return float(np.nanmean(a[floor]))
eto_m = []
for m in MONTHS:
    tmx = np.mean([prism_floor_mean("tmax", f"{y}{m:02d}") for y in range(2003, 2026)])
    tmn = np.mean([prism_floor_mean("tmin", f"{y}{m:02d}") for y in range(2003, 2026)])
    dmid = dt.date(2015, m, 15).timetuple().tm_yday
    d_ = 0.409 * math.sin(2 * math.pi * dmid / 365 - 1.39); r_ = 1 + 0.033 * math.cos(2 * math.pi * dmid / 365); p_ = math.radians(37.0)
    w_ = math.acos(-math.tan(p_) * math.tan(d_)); ra = 24 * 60 / math.pi * 0.0820 * r_ * (w_ * math.sin(p_) * math.sin(d_) + math.cos(p_) * math.cos(d_) * math.sin(w_))
    eto_m.append(0.0023 * (ra / 2.45) * ((tmx + tmn) / 2 + 17.8) * math.sqrt(max(0, tmx - tmn)))
ETo = float(np.mean(eto_m))
MN = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]
kc = pd.read_csv(os.path.join(Z, "cvhm2_landuse", "CVHM2_CropCoefficients.txt"), sep="\t").set_index("Crop")[[MN[m - 1] for m in MONTHS]].mean(1)
GROUP_CODES = {0: [1], 1: [23], 2: [3], 3: [17], 4: [4, 10, 12, 16], 5: [6, 7, 9, 11, 14, 18, 21, 22, 24], 6: [5, 13, 20], 7: [8], 8: [19], 9: [2, 15]}
tw = pd.read_csv(os.path.join(Z, "ameriflux", "tower_monthly_et.csv"))
dry_et = float(tw[(tw.site == "US-Var") & tw.month.isin(MONTHS)].et_mm.mean() / 30.4)
et_g = {g: float(np.mean([kc[c] for c in cs])) * ETo for g, cs in GROUP_CODES.items()}
et_g[2] = et_g[7] = dry_et; et_g[9] = 0.4 * ETo
print(f"Hargreaves ETo {SEASON} floor mean {ETo:.2f} mm/d; US-Var {dry_et:.2f} mm/d")
print("ET today by group (mm/d):", {GROUPS[g]: round(v, 2) for g, v in et_g.items()})

# ---- ET 1850 with the groundwater term
with rasterio.open(os.path.join(BA, f"p4_then_{P4}_mm_1km.tif")) as s:
    a = np.full((MH, MW), np.nan, np.float32)
    reproject(rasterio.band(s, 1), a, dst_transform=dst_T, dst_crs=G["crs"], resampling=Resampling.nearest, src_nodata=s.nodata, dst_nodata=np.nan)
et_p4 = np.where(a >= 0, a / NDAY, np.nan)[FJ, FI]
wells = [(f["geometry"]["coordinates"][0], f["geometry"]["coordinates"][1], 1.0 if f["properties"]["flow"] else f["properties"]["depth_m"])
         for f in fiona.open(os.path.join(BA, "p2_wsp_wells.gpkg")) if f["properties"]["on_floor"] and (f["properties"]["flow"] or (f["properties"]["depth_m"] or -1) >= 0)]
W = np.array(wells); tree = cKDTree(W[:, :2])
dist, k8 = tree.query(np.c_[X[floor], Y[floor]], k=8, distance_upper_bound=10000)
ok = np.isfinite(dist); wgt = np.where(ok, 1 / np.maximum(dist, 100) ** 2, 0)
DEPTH = np.where(ok.any(1), (wgt * np.where(ok, W[np.minimum(k8, len(W) - 1), 2], 0)).sum(1) / np.maximum(wgt.sum(1), 1e-12), np.nan)
NT = ntype[FJ, FI]
GW_EXT = {1: 3.0, 6: 3.0, 7: 3.0, 2: 6.0}                       # grassland, alkali, floodplain, oak
GW_SENS = {"low": {1: 2.0, 6: 2.0, 7: 2.0, 2: 4.0}, "high": {1: 4.0, 6: 4.0, 7: 4.0, 2: 8.0}}
KC_NATIVE = float(kc[3])
def et_1850(ext):
    e = et_p4.copy()
    if not GW:
        return e
    for t, dx in ext.items():
        m = (NT == t) & np.isfinite(DEPTH)
        e[m] = np.maximum(e[m], KC_NATIVE * ETo * np.clip(1 - DEPTH[m] / dx, 0, 1))
    return e
print(f"1850 water table: {len(W)} wells; nodes with a depth {int(np.isfinite(DEPTH).sum())} of {len(DEPTH)}; median depth {np.nanmedian(DEPTH):.1f} m; share under 3 m {np.nanmean(DEPTH < 3) * 100:.0f} %")

# ---- energy balance, vectorised over nodes
PARAM = {"water": (0.08, 0.001, 0.40), "marsh": (0.15, 0.15, 0.10), "riparian": (0.14, 1.0, 0.08), "oak": (0.17, 0.5, 0.10),
         "grass": (0.23, 0.03, 0.20), "alkali": (0.28, 0.02, 0.25), "floodplain": (0.20, 0.05, 0.15), "chaparral": (0.16, 0.4, 0.12),
         "rice": (0.15, 0.10, 0.10), "orchard": (0.17, 0.60, 0.08), "field": (0.20, 0.08, 0.12), "pasture": (0.21, 0.04, 0.12),
         "idle": (0.25, 0.02, 0.25), "urban": (0.16, 0.80, 0.20)}
NOW_CLS = {0: "water", 1: "marsh", 2: "grass", 3: "rice", 4: "orchard", 5: "field", 6: "pasture", 7: "idle", 8: "field", 9: "urban"}
THEN_CLS = {0: "grass", 1: "grass", 2: "oak", 3: "riparian", 4: "marsh", 5: "water", 6: "alkali", 7: "floodplain", 8: "chaparral"}
SIG, EPS, RHOCP, K = 5.67e-8, 0.97, 1.16e3, 0.41
ea = 0.6108 * math.exp(17.27 * Td / (Td + 237.3))
Ld = 1.24 * (ea * 10 / (Ta + 273.15)) ** (1 / 7) * SIG * (Ta + 273.15) ** 4

def params(cls):
    p = np.array([PARAM[c] if c else (np.nan,) * 3 for c in cls]); return p[:, 0], p[:, 1], p[:, 2]
def ra_stable(z0, ts, c):
    h = z0 / 0.123; d = 0.67 * h; z = np.maximum(10.0, h + 2.0) - d
    L = np.full(z0.shape, -1e6)
    for _ in range(30):
        x = np.where(L < 0, np.maximum(1e-6, 1 - 16 * z / L) ** 0.25, 1.0)
        psm = np.where(L < 0, 2 * np.log((1 + x) / 2) + np.log((1 + x * x) / 2) - 2 * np.arctan(x) + np.pi / 2, -5 * np.minimum(z / L, 1.0))
        psh = np.where(L < 0, 2 * np.log((1 + x * x) / 2), psm)
        us = K * U / np.maximum(0.05, np.log(z / z0) - psm)
        ra = (np.log(z / (0.1 * z0)) - psh) / (K * us)
        hf = RHOCP * (ts - Ta) / (c * ra)
        L = np.where(np.abs(hf) > 1e-3, -RHOCP * (Ta + 273.15) * us ** 3 / (K * 9.81 * np.where(np.abs(hf) > 1e-3, hf, 1.0)), 1e6)
    return ra
def solve(cls, et, c):
    alb, z0, gs = params(cls)
    lo, hi = np.full(len(cls), Ta - 15), np.full(len(cls), Ta + 60)
    for _ in range(40):
        ts = (lo + hi) / 2
        rn = (1 - alb) * Rs + EPS * Ld - EPS * SIG * (ts + 273.15) ** 4
        g = gs * rn; le = np.minimum(et * 28.36 * KINST, 0.95 * np.maximum(rn - g, 0))
        r = rn - g - le - RHOCP * (ts - Ta) / (c * ra_stable(np.nan_to_num(z0, nan=0.1), ts, c))
        lo = np.where(r > 0, ts, lo); hi = np.where(r > 0, hi, ts)
    out = (lo + hi) / 2
    return np.where(np.isfinite(alb) & np.isfinite(et), out, np.nan)

GR = grp[FJ, FI]
now_cls = [NOW_CLS.get(int(g)) for g in GR]; now_et = np.array([et_g.get(int(g), np.nan) for g in GR])
then_cls = [THEN_CLS.get(int(t)) if t > 0 else None for t in NT]
def skill(m_):
    k = np.isfinite(m_) & np.isfinite(OBS); dl = m_[k] - OBS[k]
    return {"bias_c": round(float(dl.mean()), 2), "rmse_c": round(float(np.sqrt((dl ** 2).mean())), 2), "r": round(float(np.corrcoef(m_[k], OBS[k])[0, 1]), 3), "n": int(k.sum())}
lo_c, hi_c = 0.05, 3.0
for _ in range(14):
    c = (lo_c + hi_c) / 2
    if skill(solve(now_cls, now_et, c))["bias_c"] > 0: hi_c = c
    else: lo_c = c
C = (lo_c + hi_c) / 2
now = solve(now_cls, now_et, C)
print(f"fitted c {C:.3f}; skill today: {skill(now)}")
k = np.isfinite(now) & np.isfinite(OBS)
bb, aa = np.polyfit(now[k], OBS[k], 1); adj = lambda v: aa + bb * v
rmse_after = float(np.sqrt(np.mean((OBS[k] - adj(now[k])) ** 2)))
print(f"adjustment observed = {aa:.2f} + {bb:.3f} model; RMSE after {rmse_after:.2f} C")

runs = {"central": GW_EXT, **GW_SENS}
res = {}
for name, ext in runs.items():
    e = et_1850(ext); th = adj(solve(then_cls, e, C))
    both = np.isfinite(th) & np.isfinite(OBS)
    res[name] = {"floor_mean_1850_c": round(float(th[both].mean()), 2), "floor_mean_today_c": round(float(OBS[both].mean()), 2),
                 "change_c": round(float((th[both] - OBS[both]).mean()), 2), "floor_et_1850_mm_d": round(float(np.nanmean(e)), 2),
                 "by_class": {TYPES[t]: [round(float(np.nanmean(e[NT == t])), 2), round(float(np.nanmean(th[(NT == t) & both])), 1)] for t in range(1, len(TYPES)) if ((NT == t) & both).sum() >= 20}}
    if name == "central":
        then_c = th
    print(name, json.dumps(res[name]))
wet_v = float(np.nanmean(OBS[GR == 1])); dry_v = float(np.nanmean(OBS[GR == 2])); wat_v = float(np.nanmean(OBS[GR == 0]))
remn = np.where(NT == 5, wat_v, np.where(np.isin(NT, [3, 4, 7]), wet_v, np.where(NT > 0, dry_v, np.nan)))
both = np.isfinite(remn) & np.isfinite(OBS) & np.isfinite(then_c)
remnant = {"floor_mean_c": round(float(remn[both].mean()), 2), "change_c": round(float((remn[both] - OBS[both]).mean()), 2), "note": "drained-state: today's native ground sits over a floor water table at a median 22 m"}
print("remnant (drained state)", remnant)

def enc(v):
    a = np.full((MH, MW), np.nan, np.float32); a[FJ, FI] = v
    return Image.fromarray(np.where(np.isfinite(a), np.clip(np.round((a - 10) * 4), 1, 255), 0).astype(np.uint8), "L")
tag = f"r3_{SEASON}" + ("" if GW else "_nogw")
enc(then_c).save(os.path.join(OUT, f"{tag}_then.png")); enc(adj(now)).save(os.path.join(OUT, f"{tag}_now.png")); enc(OBS).save(os.path.join(OUT, f"{tag}_obs.png"))
mp = os.path.join(OUT, "model_r3_meta.json")
mm = json.load(open(mp, encoding="utf-8")) if os.path.exists(mp) else {"SOURCE": __doc__.split("Writes")[0].strip(), "WHEN": "2026-11-17", "STATUS": "derived", "encoding": "C = 10 + value / 4; 0 = none", "seasons": {}}
mm["seasons"][tag] = {"days": len(DAYS), "first_days": DAYS[:3], "air": {"T_c": round(Ta, 2), "Td_c": round(Td, 2), "wind_ms": round(U, 2)}, "doy": round(float(doy)),
                      "eto_mm_d": round(ETo, 2), "et_today_by_group_mm_d": {GROUPS[g]: round(v, 2) for g, v in et_g.items()},
                      "water_table_1850": {"wells": len(W), "median_m": round(float(np.nanmedian(DEPTH)), 2), "share_under_3m": round(float(np.nanmean(DEPTH < 3)), 3)},
                      "fitted_c": round(C, 3), "skill_today": skill(now), "adjustment": {"a_c": round(float(aa), 2), "b": round(float(bb), 3), "rmse_after_c": round(rmse_after, 2)},
                      "runs": res, "remnant_drained_state": remnant}
json.dump(mm, open(mp, "w", encoding="utf-8"), indent=1)
