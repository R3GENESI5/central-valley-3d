"""Contrast adjustment of the ground-heat model (follows valley3d_heat_model.py, same day).

The energy balance with stability reproduces today's floor mean (bias 0.01 C) but its contrast between dry and watered
land is about three times the observed one (RMSE 8.78 C; dry native and idle +11 C, watered classes -5 to -7 C). Likely
causes, not separated here: 1 km MODIS pixels mix classes, and hot and cool fields exchange heat sideways (advection),
which a point energy balance does not have. The adjustment is the linear map observed = a + b * model fitted on today's
floor nodes; the same map is applied to the 1850 run, on the assumption that mixing and advection act the same way on
the 1850 land.
Writes figures/valley_3d/data/lst/model_then_adj.png (C = 10 + value / 4, 0 = none) and adds "adjustment" to
lst/model_meta.json.
WHEN 2026-11-16. STATUS derived (two fitted numbers, a and b).
Run  python src/landscape/valley3d_heat_adjust.py
"""
import json
import os

import numpy as np
from PIL import Image

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = os.path.join(REPO, "figures", "valley_3d", "data", "lst")
L = lambda f: (lambda a: np.where(a > 0, 10 + a / 4, np.nan))(np.array(Image.open(os.path.join(D, f))).astype(float))
obs, now, then = L("obs_clim.png"), L("model_now.png"), L("model_then.png")
m = np.isfinite(obs) & np.isfinite(now)
b, a = np.polyfit(now[m], obs[m], 1)
adj_now, adj_then = a + b * now, a + b * then
both = np.isfinite(adj_then) & np.isfinite(obs) & np.isfinite(now)
res = {"a_c": round(float(a), 2), "b": round(float(b), 3),
       "rmse_today_after_c": round(float(np.sqrt(np.nanmean((obs[m] - adj_now[m]) ** 2))), 2),
       "floor_mean_1850_adjusted_c": round(float(np.nanmean(adj_then[both])), 2),
       "floor_mean_today_observed_c": round(float(np.nanmean(obs[both])), 2),
       "change_1850_minus_today_c": round(float(np.nanmean(adj_then[both] - obs[both])), 2)}
print(res)
enc = np.where(np.isfinite(adj_then), np.clip(np.round((adj_then - 10) * 4), 1, 255), 0).astype(np.uint8)
Image.fromarray(enc, "L").save(os.path.join(D, "model_then_adj.png"))
mm = json.load(open(os.path.join(D, "model_meta.json"), encoding="utf-8"))
mm["adjustment"] = res | {"note": __doc__.split("Writes")[0].strip()}
json.dump(mm, open(os.path.join(D, "model_meta.json"), "w", encoding="utf-8"), indent=1)
