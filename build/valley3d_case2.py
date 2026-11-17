"""Measured layers added to the 3D page on the main California session's advice (2026-11-17, Ali: "Yes, apply all").
Station positions only; every number on the page is copied from the record named beside it, not recomputed here.
  fog      the 11 floor ISD stations of C5 (derived/m24_memory/stations.csv, grp floor)       raw/tests/2026-11-16_c5-floor-air-between-close-storms
  cold     the four hill-floor pairs of the cold-layer test (GHCN station files)               raw/tests/2026-11-17_cold-layer-hill-minus-floor
  upslope  the 45 floor and Sierra windward stations of A10b'/ERA5 bands, with their band (bands_of in
           src/landscape/snow_path_era5_bands.py)                                              raw/tests/2026-11-07_snow-path-era5-bands, 2026-10-26_mimo1-reading
Writes figures/valley_3d/data/case/case2.json.
WHEN 2026-11-17. STATUS derived (positions; numbers copied from records).
Run  python src/landscape/valley3d_case2.py
"""
import json
import os

import pandas as pd
from rasterio.warp import transform

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(REPO, "figures", "valley_3d", "data")
Z = r"Z:\Data\california_central_valley"
G = json.load(open(os.path.join(DATA, "meta.json"), encoding="utf-8"))["grid"]
def px(lat, lon):
    x, y = transform("EPSG:4326", G["crs"], [lon], [lat])
    return [round((x[0] - G["x0"]) / G["tex_m"], 1), round((G["y1"] - y[0]) / G["tex_m"], 1)]

fog = pd.read_csv(os.path.join(Z, "derived", "m24_memory", "stations.csv"), dtype={"sid": str})
fog = [[*px(r.lat, r.lon), r["name"].title()] for _, r in fog[fog.grp == "floor"].iterrows()]

ST = {"USC00040343": ("Ash Mountain", 36.4914, -118.8253, 521), "USC00040449": ("Balch Power House", 36.9092, -119.0883, 529),
      "USC00049367": ("Visalia", 36.32833, -119.2997, 92), "USC00043747": ("Hanford", 36.3158, -119.637, 72),
      "USC00045233": ("Madera", 36.9539, -120.0378, 82), "USW00093193": ("Fresno airport", 36.77999, -119.72016, 102)}
pairs = [("USC00040343", "USC00049367", -0.52, -0.98, -0.09), ("USC00040343", "USC00043747", -0.93, -1.46, -0.37),
         ("USC00040449", "USW00093193", -0.60, -0.90, -0.29), ("USC00040449", "USC00045233", -0.14, -0.44, 0.14)]
cold = [[*px(ST[h][1], ST[h][2]), ST[h][0], ST[h][3], *px(ST[f][1], ST[f][2]), ST[f][0], ST[f][3], c, lo, hi] for h, f, c, lo, hi in pairs]

def band(e, cls):
    if cls == "valley_floor": return "floor"
    return "foothills" if e < 1000 else "middle" if e < 1500 else "high"
st = pd.read_csv(os.path.join(Z, "derived", "level0", "a10bp_stations.csv"))
up = [[*px(r.lat, r.lon), r["name"].title(), round(float(r.elev)), band(r.elev, r.cls)] for _, r in st.iterrows()]

out = {"SOURCE": __doc__.split("Writes")[0].strip(), "fog_stations": fog, "cold_pairs": cold, "upslope_stations": up}
json.dump(out, open(os.path.join(DATA, "case", "case2.json"), "w", encoding="utf-8"), indent=0)
print("fog", len(fog), "| cold pairs", len(cold), "| upslope", len(up), pd.Series([u[-1] for u in up]).value_counts().to_dict())
