"""Assemble the static web map under docs/: copy rendered overlays, GeoJSON products,
weather, reports and boundaries, and write docs/data/site.json as the single manifest.

Usage: python export_site.py OUTDIR(analysis+render outputs) DOCS
"""
import os, sys, json, shutil, datetime as dt
import numpy as np
from shapely.geometry import shape, mapping
from shapely.ops import unary_union
import config as C

def simplify_fc(path, tol, props=None):
    fc = json.load(open(path))
    geoms = [shape(f["geometry"]) for f in fc["features"]]
    g = unary_union(geoms).simplify(tol, preserve_topology=True)
    return dict(type="Feature", geometry=mapping(g), properties=props or {})

def main(outdir, docs):
    data = os.path.join(docs, "data"); os.makedirs(data, exist_ok=True)
    site = dict(generated=dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%MZ"), regions={}, summary={})
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
                shutil.copy(p, os.path.join(data, f"{prod}_{region}.geojson"))
        p = os.path.join(outdir, f"summary_{region}.json")
        if os.path.exists(p):
            site["summary"][region] = json.load(open(p))
    # weather: stations as GeoJSON + forecast + recent
    w = json.load(open(os.path.join(C.SCRATCH, "weather", "weather.json")))
    st_feats = []
    for sid, s in w["stations"]["2026"].items():
        if "lat" not in s:
            continue
        prior = {y: w["stations"][y].get(sid) for y in ("2023", "2024", "2025")}
        series = {d: v for d, v in s["series"].items() if d >= "20260915"}
        st_feats.append(dict(type="Feature", geometry=dict(type="Point", coordinates=[s["lon"], s["lat"]]), properties=dict(
            id=sid, name=s["name"], elev=s["elev"], last_obs=s["last_obs"], nights_le_7=s["nights_le_7"], nights_le_2=s["nights_le_2"],
            nights_le_0=s["nights_le_0"], min_tmin=s["min_tmin"], mean_tmin=s["mean_tmin"], first_frost=s["first_frost"],
            prior={y: (None if p is None else dict(nights_le_7=p["nights_le_7"], nights_le_2=p["nights_le_2"], first_frost=p["first_frost"], mean_tmin=p["mean_tmin"])) for y, p in prior.items()},
            series=series)))
    json.dump(dict(type="FeatureCollection", features=st_feats), open(os.path.join(data, "stations.geojson"), "w"))
    json.dump(dict(gfs=w["gfs"], gfs_cycle=w["gfs_cycle"], points=w["points"], hrrr=w["hrrr"]), open(os.path.join(data, "weather.json"), "w"))
    # reports
    rep = json.load(open(os.path.join(os.path.dirname(__file__), "data", "reports.json")))
    feats = [dict(type="Feature", geometry=dict(type="Point", coordinates=[r["lon"], r["lat"]]), properties={k: v for k, v in r.items() if k not in ("lon", "lat")}) for r in rep["reports"]]
    json.dump(dict(type="FeatureCollection", retrieved=rep["retrieved"], features=feats), open(os.path.join(data, "reports.geojson"), "w"))
    # boundaries
    B = os.path.join(C.SCRATCH, "bounds")
    json.dump(dict(type="FeatureCollection", features=[simplify_fc(os.path.join(B, "king_county.geojson"), 0.002, {"name": "King County"}),
                                                       simplify_fc(os.path.join(B, "city-limits.geojson"), 0.0008, {"name": "Seattle"})]),
              open(os.path.join(data, "boundaries.geojson"), "w"))
    json.dump(site, open(os.path.join(data, "site.json"), "w"), indent=1)
    tpl = os.path.join(os.path.dirname(__file__), "templates")
    html = open(os.path.join(tpl, "index.html")).read().replace("{{LEAFLET_CSS}}", open(os.path.join(tpl, "leaflet.css")).read())
    open(os.path.join(docs, "index.html"), "w").write(html)
    tot = sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk(data) for f in fs)
    print("site manifest written; data size %.1f MB" % (tot / 1e6))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
