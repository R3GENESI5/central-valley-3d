"""Fetch CAL FIRE fire perimeters 1878-2025 for the 3D Central Valley page area, simplified to about 50 m
(Ali, 2026-11-16: "Fetch, go", after sizing: 20,802 perimeters, about 21 MB simplified).

SOURCE  CAL FIRE California_Fire_Perimeters_Time_Enabled_View FeatureServer layer 0 (services1.arcgis.com/jUJYIo9tSA7EHvfZ),
        envelope -124.5, 33.5, -117.0, 42.0 (EPSG:4326), maxAllowableOffset 0.0005 deg, fields YEAR_, FIRE_NAME, GIS_ACRES.
WRITES  Z:/Data/california_central_valley/fire_perimeters/calfire_perimeters_1878_2025_domain_gen50m.geojson (new file;
        the existing Delta-box file is untouched)
WHEN    2026-11-16. STATUS fetched (as served).
Run     py -3.12 src/landscape/valley3d_fire_fetch.py
"""
import json
import os
import time
import urllib.parse
import urllib.request

URL = "https://services1.arcgis.com/jUJYIo9tSA7EHvfZ/arcgis/rest/services/California_Fire_Perimeters_Time_Enabled_View/FeatureServer/0/query"
OUT = r"Z:\Data\california_central_valley\fire_perimeters\calfire_perimeters_1878_2025_domain_gen50m.geojson"
BASE = {"where": "1=1", "geometry": "-124.5,33.5,-117.0,42.0", "geometryType": "esriGeometryEnvelope", "inSR": "4326",
        "spatialRel": "esriSpatialRelIntersects", "outFields": "OBJECTID,YEAR_,FIRE_NAME,GIS_ACRES", "returnGeometry": "true",
        "outSR": "4326", "maxAllowableOffset": "0.0005", "orderByFields": "OBJECTID", "f": "geojson"}
assert not os.path.exists(OUT), "output exists; not overwriting"
feats, offset = [], 0
while True:
    q = dict(BASE, resultOffset=str(offset), resultRecordCount="2000")
    for attempt in range(4):
        try:
            with urllib.request.urlopen(URL + "?" + urllib.parse.urlencode(q), timeout=180) as r:
                page = json.load(r)
            break
        except Exception as e:
            print("retry", offset, e); time.sleep(10 * (attempt + 1))
    else:
        raise SystemExit("fetch failed at offset %d" % offset)
    got = page.get("features", [])
    feats += got
    print("offset", offset, "got", len(got), "total", len(feats), flush=True)
    if len(got) < 2000:
        break
    offset += 2000
    time.sleep(1)
json.dump({"type": "FeatureCollection", "features": feats,
           "source": "CAL FIRE California_Fire_Perimeters_Time_Enabled_View, fetched 2026-11-16, simplified 0.0005 deg"},
          open(OUT, "w", encoding="utf-8"))
print("wrote", OUT, os.path.getsize(OUT) // 1024, "KB,", len(feats), "perimeters")
