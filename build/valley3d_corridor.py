"""Coast-to-crest corridor for the 3D page (note raw/notes/2026-11-17_note-for-3d-page-session.md, applied on Ali's go 2026-11-17: "Yes, apply it").
Track as src/corridor/track.py defines it: start 37.80 N, 123.50 W, bearing 71 degrees, 400 km, half-width 40 km, 25 km bins;
its along/cross frame is inverted here to place a 5 km lattice of track points in page texture pixels (300 m).
Station places are the along/cross km printed in raw/tests/2026-11-17_corridor-b-c-lowest-air-and-drying (rounded there to 1 km).
Numbers on the page come from the printed output of the records named on the page, not from this file.
Writes figures/valley_3d/data/case/corridor.json.
WHEN 2026-11-17. STATUS derived (positions only).
Run  python src/landscape/valley3d_corridor.py
"""
import json
import os

import numpy as np
from rasterio.warp import transform

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(REPO, "figures", "valley_3d", "data")
G = json.load(open(os.path.join(DATA, "meta.json"), encoding="utf-8"))["grid"]
LAT0, LON0, BEARING, LENGTH_KM, HALF_KM, BIN_KM, KM = 37.80, -123.50, 71.0, 400.0, 40.0, 25.0, 111.195   # as track.py


def latlon(a, c):
    """Inverse of track.along_cross: along a km, cross c km (positive to the right of the bearing, the south side, as track.py signs it)."""
    b = np.radians(BEARING)
    dx, dy = a * np.sin(b) + c * np.cos(b), a * np.cos(b) - c * np.sin(b)
    return LAT0 + dy / KM, LON0 + dx / (KM * np.cos(np.radians(LAT0)))


def px(a, c):
    lat, lon = latlon(np.asarray(a, float), np.asarray(c, float))
    x, y = transform("EPSG:4326", G["crs"], list(np.ravel(lon)), list(np.ravel(lat)))
    return [[round((xi - G["x0"]) / G["tex_m"], 1), round((G["y1"] - yi) / G["tex_m"], 1)] for xi, yi in zip(x, y)]


A = np.arange(0, LENGTH_KM + 0.1, 5.0); C = np.arange(-HALF_KM, HALF_KM + 0.1, 5.0)
aa, cc = np.meshgrid(A, C, indexing="ij")
lattice = px(aa.ravel(), cc.ravel())
# along km, cross km, name, group: as printed in corridor-b-c (stations) and its buoy line
ST = [(53, 24, "Buoy 46026", "sea"), (88, 51, "San Francisco Intl", "rim"), (102, 45, "Oakland Intl", "rim"),
      (147, -4, "Travis AFB", "gap"), (191, -40, "Sacramento Intl", "floor"), (192, -17, "Sacramento Executive", "floor"),
      (206, -31, "McClellan", "floor"), (192, 56, "Stockton", "floor"), (285, -75, "Blue Canyon, 1,611 m", "slope"),
      (372, -73, "Reno, 1,343 m", "east")]
stations = [p + [n, g, a, c] for (a, c, n, g), p in zip(ST, px([s[0] for s in ST], [s[1] for s in ST]))]
out = {"SOURCE": __doc__.split("Writes")[0].strip(), "step_km": 5.0, "na": len(A), "nc": len(C), "bin_km": BIN_KM,
       "lattice": lattice, "stations": stations}
json.dump(out, open(os.path.join(DATA, "case", "corridor.json"), "w", encoding="utf-8"), separators=(",", ":"))
inside = sum(0 <= x < G["tw"] and 0 <= y < G["th"] for x, y in lattice)
print("lattice", len(A), "x", len(C), "| points inside the page grid", inside, "of", len(lattice))
for s in stations: print(s)
