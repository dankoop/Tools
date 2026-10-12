"""Read Sentinel-2 L2A scenes (cloud-optimized GeoTIFFs on AWS) into analysis-ready
per-scene band stacks on a fixed UTM grid, with a quality mask derived from the scene
classification layer (SCL) plus spectral cloud/shadow tests.

Output per scene: <out>/<date>_<id>.npz with
    bands  uint16 [nband, rows, cols]  raw L2A digital numbers (0 = nodata)
    scl    uint8  [rows, cols]
    valid  bool   [rows, cols]          clear-sky land observation
Multiple tiles acquired the same day are mosaicked (first valid pixel wins, clearest tile first).

Usage: python ingest_s2.py REGION SCENES.json OUTDIR [--min-date YYYY-MM-DD] [--res-override M]
"""
import os, sys, json, argparse, math, concurrent.futures as cf
import numpy as np, rasterio
from rasterio.windows import Window
from rasterio.enums import Resampling
from pyproj import Transformer
from scipy import ndimage
import config as C

os.environ.update(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
                  GDAL_HTTP_MULTIRANGE="YES", GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES",
                  VSI_CACHE="TRUE", VSI_CACHE_SIZE="268435456", GDAL_HTTP_MAX_RETRY="5", GDAL_HTTP_RETRY_DELAY="2")

def region_grid(name, res=None):
    """Return (transform, width, height) for a region on a UTM grid snapped to `res`."""
    r = C.REGIONS[name]; res = res or r["res"]
    tr = Transformer.from_crs(4326, C.UTM, always_xy=True)
    lon0, lat0, lon1, lat1 = r["bbox"]
    xs, ys = zip(*[tr.transform(x, y) for x, y in [(lon0, lat0), (lon0, lat1), (lon1, lat0), (lon1, lat1)]])
    x0 = math.floor(min(xs) / res) * res; x1 = math.ceil(max(xs) / res) * res
    y0 = math.floor(min(ys) / res) * res; y1 = math.ceil(max(ys) / res) * res
    transform = rasterio.transform.from_origin(x0, y1, res, res)
    return transform, int((x1 - x0) / res), int((y1 - y0) / res)

def read_band_to_grid(href, transform, width, height, res, resampling):
    """Read a COG into the target grid. Returns uint16 array (0 where outside the tile)."""
    out = np.zeros((height, width), np.uint16)
    with rasterio.open("/vsicurl/" + href) as ds:
        b = ds.bounds
        gx0, gy1 = transform.c, transform.f
        gx1, gy0 = gx0 + width * res, gy1 - height * res
        ix0, ix1 = max(gx0, b.left), min(gx1, b.right)
        iy0, iy1 = max(gy0, b.bottom), min(gy1, b.top)
        if ix1 <= ix0 or iy1 <= iy0:
            return out
        # target-grid indices of the intersection
        c0 = int(round((ix0 - gx0) / res)); c1 = int(round((ix1 - gx0) / res))
        r0 = int(round((gy1 - iy1) / res)); r1 = int(round((gy1 - iy0) / res))
        src_res = ds.res[0]
        sc0 = (ix0 - b.left) / src_res; sr0 = (b.top - iy1) / src_res
        win = Window(sc0, sr0, (ix1 - ix0) / src_res, (iy1 - iy0) / src_res)
        data = ds.read(1, window=win, out_shape=(r1 - r0, c1 - c0), resampling=resampling, boundless=True, fill_value=0)
        out[r0:r1, c0:c1] = data
    return out

def scene_quality_mask(bands, scl):
    """Clear-sky land mask: SCL keep-classes, spectral cloud/shadow tests, 2-px cloud dilation."""
    return quality_mask_from(bands[C.BAND_ORDER.index("blue")], bands[C.BAND_ORDER.index("nir")],
                             bands[C.BAND_ORDER.index("swir16")], (bands == 0).any(axis=0), scl)

def quality_mask_from(blue_dn, nir_dn, swir_dn, nodata, scl):
    """Same mask from the three bands it needs (lets callers skip building the full band stack)."""
    blue = (blue_dn.astype(np.float32) - C.BOA_OFFSET) / C.BOA_SCALE
    nir = (nir_dn.astype(np.float32) - C.BOA_OFFSET) / C.BOA_SCALE
    swir = (swir_dn.astype(np.float32) - C.BOA_OFFSET) / C.BOA_SCALE
    nodata = nodata | (scl == 0)
    cloudy = np.isin(scl, [1, 2, 3, 8, 9, 10, 11]) | (blue > 0.22) | ((blue > 0.12) & (swir > 0.35))
    cloudy = ndimage.binary_dilation(cloudy, iterations=2)
    shadow = (nir < 0.06) & (scl != C.SCL_WATER)
    keep = np.isin(scl, list(C.SCL_KEEP))
    return keep & ~cloudy & ~shadow & ~nodata

def ingest_scene_group(items, transform, width, height, res, out_path):
    """Mosaic all tiles for one date into the grid and save."""
    if os.path.exists(out_path):
        return out_path + " (cached)"
    items = sorted(items, key=lambda i: i["cloud"])      # clearest tile first
    bands = np.zeros((len(C.BAND_ORDER), height, width), np.uint16)
    scl = np.zeros((height, width), np.uint8)
    valid = np.zeros((height, width), bool)
    for it in items:
        jobs = {b: (it["assets"][b], Resampling.average if C.BANDS[b] < res else Resampling.bilinear)
                for b in C.BAND_ORDER}
        with cf.ThreadPoolExecutor(6) as ex:
            futs = {b: ex.submit(read_band_to_grid, h, transform, width, height, res, rs) for b, (h, rs) in jobs.items()}
            fscl = ex.submit(read_band_to_grid, it["assets"]["scl"], transform, width, height, res, Resampling.nearest)
            tb = np.stack([futs[b].result() for b in C.BAND_ORDER])
            ts = fscl.result().astype(np.uint8)
        tv = scene_quality_mask(tb, ts)
        fill = ~valid & (tb != 0).all(axis=0) & (bands == 0).any(axis=0)   # pixels not yet covered
        take = (tv & ~valid) | fill
        bands[:, take] = tb[:, take]; scl[take] = ts[take]; valid[take] = tv[take] | valid[take]
        valid |= tv & (bands != 0).all(axis=0)
    np.savez_compressed(out_path, bands=bands, scl=scl, valid=valid,
                        ids=np.array([i["id"] for i in items]), band_order=np.array(C.BAND_ORDER))
    return f"{out_path} valid={valid.mean():.3f}"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("region"); ap.add_argument("scenes"); ap.add_argument("outdir")
    ap.add_argument("--min-date", default="2026-09-01"); ap.add_argument("--max-date", default="2099-01-01")
    ap.add_argument("--res", type=int, default=None, help="override grid resolution (m)")
    ap.add_argument("--max-cloud", type=float, default=95)
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    transform, width, height = region_grid(a.region, a.res)
    res = a.res or C.REGIONS[a.region]["res"]
    json.dump({"transform": list(transform)[:6], "width": width, "height": height, "res": res, "crs": C.UTM},
              open(os.path.join(a.outdir, "grid.json"), "w"))
    scenes = [s for s in json.load(open(a.scenes)) if a.min_date <= s["date"] <= a.max_date
              and s["cloud"] < a.max_cloud and not (s["nodata"] > 90 and s["cloud"] > 60)]
    groups = {}
    for s in scenes:
        groups.setdefault(s["date"], []).append(s)
    print(f"{a.region}: grid {width}x{height} @ {res} m; {len(groups)} dates, {len(scenes)} scenes")
    jobs = [(items, transform, width, height, res, os.path.join(a.outdir, f"{d}.npz")) for d, items in sorted(groups.items())]
    with cf.ThreadPoolExecutor(a.workers) as ex:
        for r in ex.map(lambda j: ingest_scene_group(*j), jobs):
            print(r, flush=True)

if __name__ == "__main__":
    main()
