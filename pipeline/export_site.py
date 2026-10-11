"""Assemble the static web map under docs/: copy rendered overlays, GeoJSON products,
weather, reports and boundaries, and write docs/data/site.json as the single manifest.

Usage: python export_site.py OUTDIR(analysis+render outputs) DOCS
"""
import os, sys, json, shutil, datetime as dt
import numpy as np
from shapely.geometry import shape, mapping
from shapely.ops import unary_union
import config as C

def clean(o):
    """Recursively replace NaN/inf floats with None so the JSON is valid for browsers."""
    import math
    if isinstance(o, float):
        return None if (math.isnan(o) or math.isinf(o)) else o
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, list):
        return [clean(v) for v in o]
    return o

def dump(o, path):
    json.dump(clean(o), open(path, "w"), allow_nan=False)

def glob_files(d, prefixes):
    return [os.path.join(d, f) for f in os.listdir(d) if f.startswith(prefixes)]

def simplify_fc(path, tol, props=None):
    fc = json.load(open(path))
    geoms = [shape(f["geometry"]) for f in fc["features"]]
    g = unary_union(geoms).simplify(tol, preserve_topology=True)
    return dict(type="Feature", geometry=mapping(g), properties=props or {})

def main(outdir, docs):
    data = os.path.join(docs, "data"); os.makedirs(data, exist_ok=True)
    site = dict(generated=dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%MZ"), regions={}, summary={})
    merged = {}
    for region in ("seattle", "county"):
        src = os.path.join(outdir, "render", region)
        if not os.path.exists(os.path.join(src, "layers.json")):
            continue
        dst = os.path.join(data, region); os.makedirs(dst, exist_ok=True)
        for f in os.listdir(src):
            shutil.copy(os.path.join(src, f), dst)
        site["regions"][region] = json.load(open(os.path.join(src, "layers.json")))
        for prod in ("clusters", "destinations"):
            p = os.path.join(outdir, f"{prod}_{region}.geojson")
            if os.path.exists(p):
                merged.setdefault(prod, []).extend(json.load(open(p))["features"])
        p = os.path.join(outdir, f"summary_{region}.json")
        if os.path.exists(p):
            site["summary"][region] = clean(json.load(open(p)))
    # one ranking across both regions
    for prod, out_name in (("clusters", "hotspots.geojson"), ("destinations", "destinations.geojson")):
        feats = sorted(merged.get(prod, []), key=lambda f: -f["properties"]["visit_score"])
        for r, f in enumerate(feats, 1):
            f["properties"]["rank"] = r
        dump(dict(type="FeatureCollection", features=feats), os.path.join(data, out_name))
    for stale in glob_files(data, ("clusters_", "destinations_")):
        os.remove(stale)
    # weather: stations as GeoJSON + forecast + recent
    w = json.load(open(os.path.join(C.SCRATCH, "weather", "weather.json")))
    st_feats = []
    def window_stats(st, year, md_end):
        """Cold-night counts for Sep 15 .. md_end of a prior year, so the comparison matches this year's window."""
        if st is None:
            return None
        tm = {d: v[0] for d, v in st["series"].items() if v[0] is not None and "0915" <= d[4:] <= md_end}
        if not tm:
            return None
        return dict(nights_le_7=sum(v <= 7 for v in tm.values()), nights_le_2=sum(v <= 2 for v in tm.values()),
                    mean_tmin=round(sum(tm.values()) / len(tm), 1), first_frost=st.get("first_frost"), n=len(tm))
    for sid, s in w["stations"]["2026"].items():
        if "lat" not in s:
            continue
        prior = {y: window_stats(w["stations"][y].get(sid), y, s["last_obs"][4:]) for y in ("2023", "2024", "2025")}
        series = {d: v for d, v in s["series"].items() if d >= "20260915"}
        st_feats.append(dict(type="Feature", geometry=dict(type="Point", coordinates=[s["lon"], s["lat"]]), properties=dict(
            id=sid, name=s["name"], elev=s["elev"], last_obs=s["last_obs"], nights_le_7=s["nights_le_7"], nights_le_2=s["nights_le_2"],
            nights_le_0=s["nights_le_0"], min_tmin=s["min_tmin"], mean_tmin=s["mean_tmin"], first_frost=s["first_frost"],
            prior=prior,
            series=series)))
    dump(dict(type="FeatureCollection", features=st_feats), os.path.join(data, "stations.geojson"))
    dump(dict(gfs=w["gfs"], gfs_cycle=w["gfs_cycle"], points=w["points"], hrrr=w["hrrr"]), os.path.join(data, "weather.json"))
    shutil.copy(os.path.join(os.path.dirname(__file__), "data", "places.json"), os.path.join(data, "places.json"))
    # sharp 10 m base tiles and Overture context (roads, water)
    bsrc = os.path.join(outdir, "base")
    if os.path.exists(os.path.join(bsrc, "base.json")):
        bdst = os.path.join(data, "base"); os.makedirs(bdst, exist_ok=True)
        for f in os.listdir(bsrc):
            shutil.copy(os.path.join(bsrc, f), bdst)
        site["base"] = json.load(open(os.path.join(bsrc, "base.json")))
    csrc = os.path.join(C.SCRATCH, "context")
    if os.path.exists(csrc):
        cdst = os.path.join(data, "context"); os.makedirs(cdst, exist_ok=True)
        for f in ("roads_major.geojson", "roads_minor.geojson", "water.geojson"):
            if os.path.exists(os.path.join(csrc, f)):
                shutil.copy(os.path.join(csrc, f), cdst)
        site["context"] = sorted(os.listdir(cdst))
    # reports
    rep = json.load(open(os.path.join(os.path.dirname(__file__), "data", "reports.json")))
    feats = [dict(type="Feature", geometry=dict(type="Point", coordinates=[r["lon"], r["lat"]]), properties={k: v for k, v in r.items() if k not in ("lon", "lat")}) for r in rep["reports"]]
    dump(dict(type="FeatureCollection", retrieved=rep["retrieved"], features=feats), os.path.join(data, "reports.geojson"))
    # boundaries
    B = os.path.join(C.SCRATCH, "bounds")
    json.dump(dict(type="FeatureCollection", features=[simplify_fc(os.path.join(B, "king_county.geojson"), 0.002, {"name": "King County"}),
                                                       simplify_fc(os.path.join(B, "city-limits.geojson"), 0.0008, {"name": "Seattle"})]),
              open(os.path.join(data, "boundaries.geojson"), "w"))
    dump(site, os.path.join(data, "site.json"))
    tpl = os.path.join(os.path.dirname(__file__), "templates")
    html = open(os.path.join(tpl, "index.html")).read().replace("{{LEAFLET_CSS}}", open(os.path.join(tpl, "leaflet.css")).read()).replace("{{LEAFLET_JS}}", open(os.path.join(tpl, "leaflet.js")).read())
    open(os.path.join(docs, "index.html"), "w").write(html)
    tot = sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk(data) for f in fs)
    print("site manifest written; data size %.1f MB" % (tot / 1e6))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
