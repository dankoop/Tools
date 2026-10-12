"""Preset regions with summaries: Seattle neighbourhoods (from the city boundaries file,
grouped by district) and hand-defined King County areas. Each region gets stage shares,
colour score, predicted peak, weekly series and forecast from the analysis grids.

Usage: python regions.py OUT.geojson
"""
import os, sys, json, datetime as dt, collections
import numpy as np
from rasterio import features
from rasterio.transform import Affine
from shapely.geometry import shape, mapping, LineString, Polygon, box
from shapely.ops import unary_union
import shapely
from pyproj import Transformer
import config as C
from analyze import STAGE_NAMES

CACHE = C.SCRATCH
TODAY = dt.date.today()

# King County areas: corridors as buffered lines, districts as polygons (lon, lat)
COUNTY_AREAS = {
    "Eastside (Bellevue, Kirkland, Redmond)": Polygon([(-122.25, 47.52), (-122.25, 47.76), (-122.03, 47.76), (-122.03, 47.52)]),
    "Issaquah Alps (Cougar, Squak, Tiger)": Polygon([(-122.16, 47.46), (-122.16, 47.56), (-121.92, 47.56), (-121.92, 47.46)]),
    "Snoqualmie Valley (Duvall to North Bend)": LineString([(-121.985, 47.745), (-121.915, 47.65), (-121.89, 47.57), (-121.83, 47.53), (-121.78, 47.495)]).buffer(0.035),
    "I-90 corridor (North Bend to the Pass)": LineString([(-121.78, 47.49), (-121.66, 47.45), (-121.55, 47.42), (-121.45, 47.41), (-121.40, 47.42)]).buffer(0.03),
    "Highway 2 (Skykomish valley)": LineString([(-121.70, 47.72), (-121.55, 47.74), (-121.45, 47.72), (-121.36, 47.71), (-121.27, 47.73)]).buffer(0.03),
    "Cedar River and Maple Valley": LineString([(-122.20, 47.48), (-122.10, 47.42), (-122.02, 47.39), (-121.95, 47.37)]).buffer(0.03),
    "Green River and Enumclaw plateau": Polygon([(-122.12, 47.16), (-122.12, 47.33), (-121.88, 47.33), (-121.88, 47.16)]),
    "Greenwater and Highway 410": LineString([(-121.75, 47.25), (-121.68, 47.16), (-121.64, 47.10), (-121.60, 47.03)]).buffer(0.03),
    "South King (Renton to Federal Way)": Polygon([(-122.40, 47.28), (-122.40, 47.52), (-122.15, 47.52), (-122.15, 47.28)]),
    "North lake (Shoreline, Kenmore, Bothell, Woodinville)": Polygon([(-122.40, 47.73), (-122.40, 47.79), (-122.10, 47.79), (-122.10, 47.73)]),
    "Vashon and Maury Islands": Polygon([(-122.55, 47.33), (-122.55, 47.53), (-122.38, 47.53), (-122.38, 47.33)]),
}

def load_region(name):
    z = np.load(os.path.join(CACHE, "out", f"analysis_{name}.npz"))
    T = Affine(*z["transform"][:6]); return z, T

def summarize(z, mask):
    veg = z["veg"] & mask; decid = (z["decid"].astype(np.float32) - 1) / 254
    vv = veg & (decid >= 0.25)
    n = int(vv.sum())
    if n < 50:
        return None
    q = lambda a, lo, hi: np.where(a > 0, (a.astype(np.float32) - 1) / 254 * (hi - lo) + lo, np.nan)
    stage = z["stage"]; timing = z["timing"]; cs = q(z["cs_now"], 0, 1); d2p = q(z["days_to_peak"], -21, 35)
    res = int(z["res"]); out = dict(
        tree_ha=round(n * res * res / 1e4, 1),
        stage_fraction={STAGE_NAMES[i]: round(float((stage[vv] == i).mean()), 3) for i in range(1, 5)},
        intensity=round(float(np.nanmean(cs[vv])), 3),
        timing_fraction={k: round(float((timing[vv] == i).mean()), 3) for i, k in enumerate(["past peak", "peak now", "next 2 weeks", "2-4 weeks", "later"], 1)},
        series={str(l): round(float(np.nanmean(q(z[f"wk|{l}|cs"], 0, 1)[vv])), 3) for l in z["bins"]},
        forecast={}, forecast_stage={})
    dv = d2p[vv]; dv = dv[np.isfinite(dv)]
    if dv.size:
        med = float(np.median(dv)); pf = out["stage_fraction"]["peak"]; lf = out["stage_fraction"]["past peak"]
        if pf >= 0.35: med = min(med, 0.0)
        if lf >= 0.3: med = min(med, -14.0)
        out["days_to_peak"] = round(med, 1); out["peak_date"] = str(TODAY + dt.timedelta(days=int(round(med))))
        out["peak_window_days"] = int(round(5 + 12 * (1 - float(np.nanmean((z["confidence"].astype(np.float32) - 1) / 254)))))
    for rec in z["fc_bins"]:
        lab = str(rec).split("|")[0]
        out["forecast"][lab] = round(float(np.nanmean(q(z[f"fc|{lab}|cs"], 0, 1)[vv])), 3)
        st = z[f"fc|{lab}|stage"]; out["forecast_stage"][lab] = {STAGE_NAMES[i]: round(float((st[vv] == i).mean()), 3) for i in range(1, 5)}
    return out

def main(out_path):
    tr = Transformer.from_crs(4326, C.UTM, always_xy=True); to_ll = Transformer.from_crs(C.UTM, 4326, always_xy=True)
    utm = lambda g: shapely.transform(g, lambda xy: np.column_stack(tr.transform(xy[:, 0], xy[:, 1])))
    ll = lambda g: shapely.transform(g, lambda xy: np.column_stack(to_ll.transform(xy[:, 0], xy[:, 1])))
    feats = []
    # --- Seattle neighbourhoods grouped by district
    z, T = load_region("seattle"); H, W = z["veg"].shape
    nh = json.load(open(os.path.join(CACHE, "bounds", "neighborhoods.geojson")))
    groups = collections.defaultdict(list)
    for f in nh["features"]:
        p = f["properties"]
        if p.get("city") != "Seattle": continue
        groups[p.get("nhood") or p.get("name")].append(shape(f["geometry"]))
    for name, geoms in sorted(groups.items()):
        g = unary_union(geoms).buffer(0); gu = utm(g)
        m = features.rasterize([(mapping(gu), 1)], out_shape=(H, W), transform=T, dtype=np.uint8).astype(bool)
        s = summarize(z, m)
        if s is None: continue
        c = g.centroid
        feats.append(dict(type="Feature", geometry=mapping(g.simplify(0.0003)), properties=dict(id="sea-" + name.lower().replace(" ", "-").replace("/", "-"), name=name, group="Seattle neighbourhood", lon=round(c.x, 5), lat=round(c.y, 5), bbox=[round(v, 5) for v in g.bounds], **s)))
    print(len(feats), "Seattle neighbourhoods summarised")
    # --- county areas (clipped to the county, water excluded by the veg mask)
    z, T = load_region("county"); H, W = z["veg"].shape
    kc = shape(json.load(open(os.path.join(CACHE, "bounds", "king_county.geojson")))["features"][0]["geometry"])
    for name, g in COUNTY_AREAS.items():
        g = g.intersection(kc)
        m = features.rasterize([(mapping(utm(g)), 1)], out_shape=(H, W), transform=T, dtype=np.uint8).astype(bool)
        s = summarize(z, m)
        if s is None: print("skip", name); continue
        c = g.centroid
        feats.append(dict(type="Feature", geometry=mapping(g.simplify(0.001)), properties=dict(id="kc-" + name.split(" (")[0].lower().replace(" ", "-"), name=name, group="King County area", lon=round(c.x, 5), lat=round(c.y, 5), bbox=[round(v, 5) for v in g.bounds], **s)))
    json.dump(dict(type="FeatureCollection", features=feats), open(out_path, "w"))
    print(len(feats), "regions ->", out_path)
    for f in feats[-11:]:
        p = f["properties"]; print(f"  {p['name'][:40]:40s} trees {p['tree_ha']:>7} ha  peak {p['stage_fraction']['peak']:.2f} turning {p['stage_fraction']['turning']:.2f}  peak_date {p.get('peak_date')}")

if __name__ == "__main__":
    main(sys.argv[1])
