"""Build an inventory of Sentinel-2 L2A scenes (Element84 sentinel-cogs bucket) for the
MGRS tiles covering King County, with per-scene STAC metadata (cloud cover, nodata, etc.).

Usage: python inventory.py OUT.json YEAR MONTHS... (e.g. 2026 9 10)
"""
import json, re, sys, concurrent.futures as cf, urllib.request

BUCKET = "https://sentinel-cogs.s3.us-west-2.amazonaws.com"
TILES = ["10/T/ET", "10/T/FT"]

def list_prefixes(prefix):
    url = f"{BUCKET}/?list-type=2&prefix={prefix}&delimiter=/"
    xml = urllib.request.urlopen(url, timeout=60).read().decode()
    return re.findall(r"<Prefix>([^<]+)</Prefix>", xml)

def fetch_item(scene_prefix):
    sid = scene_prefix.rstrip("/").split("/")[-1]
    url = f"{BUCKET}/{scene_prefix}{sid}.json"
    try:
        d = json.load(urllib.request.urlopen(url, timeout=60))
    except Exception as e:
        return {"id": sid, "error": str(e)}
    p = d["properties"]
    return {
        "id": sid, "tile": sid.split("_")[1], "datetime": p.get("datetime"),
        "date": p.get("datetime", "")[:10], "cloud": p.get("eo:cloud_cover"),
        "nodata": p.get("s2:nodata_pixel_percentage"), "veg": p.get("s2:vegetation_percentage"),
        "snow": p.get("s2:snow_ice_percentage"), "baseline": p.get("s2:processing_baseline"),
        "bbox": d.get("bbox"), "prefix": scene_prefix,
        "assets": {k: v["href"] for k, v in d["assets"].items() if k in
                   ("blue", "green", "red", "rededge1", "rededge2", "rededge3", "nir", "nir08",
                    "swir16", "swir22", "scl", "visual")},
    }

def main(out, year, months):
    scenes = []
    for t in TILES:
        for m in months:
            for p in list_prefixes(f"sentinel-s2-l2a-cogs/{t}/{year}/{m}/"):
                if re.search(r"/S2[ABC]_\w+_L2A/$", p):
                    scenes.append(p)
    with cf.ThreadPoolExecutor(16) as ex:
        items = list(ex.map(fetch_item, scenes))
    items = [i for i in items if "error" not in i]
    items.sort(key=lambda i: (i["datetime"], i["tile"]))
    json.dump(items, open(out, "w"), indent=1)
    print(f"{len(items)} scenes -> {out}")
    for i in items:
        print(f'{i["id"]:32s} cloud={i["cloud"]:5.1f} nodata={i["nodata"]:5.1f} veg={i["veg"]:5.1f}')

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), [int(m) for m in sys.argv[3:]])
