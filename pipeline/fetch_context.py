"""Fetch orientation context from Overture Maps (public parquet on AWS, no key):
roads and water polygons for King County, written as simplified GeoJSON.

Row groups are pruned with the parquet bbox statistics, so only a small slice of the
global files is downloaded.

Usage: python fetch_context.py OUTDIR [--release 2026-09-23.1]
"""
import os, sys, json, argparse, time, concurrent.futures as cf
import pyarrow.parquet as pq, pyarrow.compute as pc, pyarrow as pa, s3fs
from shapely import wkb, box
from shapely.geometry import mapping
from shapely.ops import unary_union
import shapely

BBOX = (-122.56, 47.06, -121.04, 47.80)          # King County with a small margin
SEATTLE_BBOX = (-122.46, 47.47, -122.20, 47.76)  # city + 2 km for minor roads
MAJOR = {"motorway", "trunk", "primary", "secondary", "tertiary"}
MINOR = {"residential", "unclassified", "living_street"}
WATER = {"lake", "reservoir", "river", "ocean", "canal", "bay", "sea", "pond", "lagoon", "stream"}

def rg_hits(md, idx, bbox):
    out = []
    for rg in range(md.num_row_groups):
        r = md.row_group(rg)
        s = {k: r.column(idx[k]).statistics for k in ("bbox.xmin", "bbox.xmax", "bbox.ymin", "bbox.ymax")}
        if any(v is None for v in s.values()):
            out.append(rg); continue
        if s["bbox.xmin"].max > bbox[0] and s["bbox.xmax"].min < bbox[2] and s["bbox.ymin"].max > bbox[1] and s["bbox.ymax"].min < bbox[3]:
            out.append(rg)
    return out

def scan(fs, path, columns, bbox, flt):
    pf = pq.ParquetFile(fs.open(path, "rb"))
    md = pf.metadata
    idx = {md.row_group(0).column(i).path_in_schema: i for i in range(md.row_group(0).num_columns)}
    hits = rg_hits(md, idx, bbox)
    if not hits:
        return None
    t = pf.read_row_groups(hits, columns=columns)
    m = (pc.field("bbox", "xmax") > bbox[0]) & (pc.field("bbox", "xmin") < bbox[2]) & (pc.field("bbox", "ymax") > bbox[1]) & (pc.field("bbox", "ymin") < bbox[3])
    t = t.filter(m)
    if flt is not None:
        t = t.filter(flt)
    return t if t.num_rows else None

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("outdir"); ap.add_argument("--release", default="2026-09-23.1")
    a = ap.parse_args(); os.makedirs(a.outdir, exist_ok=True)
    fs = s3fs.S3FileSystem(anon=True, client_kwargs={"region_name": "us-west-2"})
    base = f"overturemaps-us-west-2/release/{a.release}"
    clip = box(*BBOX)

    # ---------------- roads
    t0 = time.time()
    files = fs.ls(f"{base}/theme=transportation/type=segment/")
    flt = (pc.field("subtype") == "road") & pc.field("class").isin(list(MAJOR | MINOR))
    with cf.ThreadPoolExecutor(8) as ex:
        tabs = [t for t in ex.map(lambda f: scan(fs, f, ["id", "subtype", "class", "names", "geometry", "bbox"], BBOX, flt), files) if t is not None]
    tab = pa.concat_tables(tabs)
    print(f"roads: {tab.num_rows} segments from {len(tabs)} files in {time.time()-t0:.0f}s", flush=True)
    major, minor = [], []
    sea = box(*SEATTLE_BBOX)
    for cls, nm, g in zip(tab["class"].to_pylist(), tab["names"].to_pylist(), tab["geometry"].to_pylist()):
        geom = wkb.loads(g)
        if not geom.intersects(clip):
            continue
        geom = geom.intersection(clip)
        name = (nm or {}).get("primary") if isinstance(nm, dict) else None
        if cls in MAJOR:
            major.append(dict(type="Feature", geometry=mapping(shapely.set_precision(geom.simplify(0.00008), 0.00001)), properties=dict(c=cls[:3], n=name if cls in ("motorway", "trunk", "primary") else None)))
        elif geom.intersects(sea):
            minor.append(dict(type="Feature", geometry=mapping(shapely.set_precision(geom.intersection(sea).simplify(0.00005), 0.00001)), properties=dict(c=cls[:3])))
    for name, feats in (("roads_major.geojson", major), ("roads_minor.geojson", minor)):
        p = os.path.join(a.outdir, name)
        json.dump(dict(type="FeatureCollection", features=feats), open(p, "w"), separators=(",", ":"))
        print(name, len(feats), "features,", os.path.getsize(p) // 1024, "KB", flush=True)

    # ---------------- water
    t0 = time.time()
    files = fs.ls(f"{base}/theme=base/type=water/")
    flt = pc.field("subtype").isin(list(WATER))
    with cf.ThreadPoolExecutor(8) as ex:
        tabs = [t for t in ex.map(lambda f: scan(fs, f, ["id", "subtype", "class", "names", "geometry", "bbox"], BBOX, flt), files) if t is not None]
    tab = pa.concat_tables(tabs)
    print(f"water: {tab.num_rows} features from {len(tabs)} files in {time.time()-t0:.0f}s", flush=True)
    feats = []
    for st, nm, g in zip(tab["subtype"].to_pylist(), tab["names"].to_pylist(), tab["geometry"].to_pylist()):
        geom = wkb.loads(g)
        if geom.geom_type not in ("Polygon", "MultiPolygon") or not geom.intersects(clip):
            continue
        geom = geom.intersection(clip)
        if geom.is_empty or geom.area < 2e-6:          # skip tiny ponds (~2 ha)
            continue
        name = (nm or {}).get("primary") if isinstance(nm, dict) else None
        feats.append(dict(type="Feature", geometry=mapping(shapely.set_precision(geom.simplify(0.00006), 0.00001)), properties=dict(s=st, n=name if geom.area > 2e-4 else None)))
    p = os.path.join(a.outdir, "water.geojson")
    json.dump(dict(type="FeatureCollection", features=feats), open(p, "w"), separators=(",", ":"))
    print("water.geojson", len(feats), "features,", os.path.getsize(p) // 1024, "KB", flush=True)

if __name__ == "__main__":
    main()
