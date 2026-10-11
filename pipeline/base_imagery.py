"""Build a sharp 10 m Sentinel-2 true-colour base for the whole county from the clearest
recent passes (latest clear pixel wins), with a gentle contrast stretch and unsharp mask,
reprojected to Web Mercator and cut into <= 2048 px WebP tiles for the map.

Usage: python base_imagery.py SCENES.json OUTDIR [--min-date 2026-09-15] [--max-cloud 60]
Writes OUTDIR/base/tile_r_c.webp and OUTDIR/base/base.json (bounds per tile).
"""
import os, sys, json, argparse, math, concurrent.futures as cf
import numpy as np, rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject, calculate_default_transform
from rasterio.transform import Affine
from scipy import ndimage
from PIL import Image
import config as C
from ingest_s2 import region_grid, read_band_to_grid, scene_quality_mask

os.environ.update(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif", GDAL_HTTP_MULTIRANGE="YES",
                  GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES", VSI_CACHE="TRUE", GDAL_HTTP_MAX_RETRY="5", GDAL_HTTP_RETRY_DELAY="2")

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("scenes"); ap.add_argument("outdir")
    ap.add_argument("--min-date", default="2026-09-15"); ap.add_argument("--max-cloud", type=float, default=60); ap.add_argument("--tile", type=int, default=2048)
    a = ap.parse_args(); out = os.path.join(a.outdir, "base"); os.makedirs(out, exist_ok=True)
    T, W, H = region_grid("county", 10); res = 10
    scenes = [s for s in json.load(open(a.scenes)) if s["date"] >= a.min_date and s["cloud"] <= a.max_cloud and s["nodata"] < 90]
    # latest first so the newest clear pixel wins
    dates = sorted({s["date"] for s in scenes}, reverse=True)
    rgb = np.zeros((3, H, W), np.uint16); filled = np.zeros((H, W), bool)
    for d in dates:
        for s in sorted([x for x in scenes if x["date"] == d], key=lambda x: x["cloud"]):
            with cf.ThreadPoolExecutor(5) as ex:
                fb = {b: ex.submit(read_band_to_grid, s["assets"][b], T, W, H, res, Resampling.nearest) for b in ("red", "green", "blue", "nir", "swir16")}
                fs = ex.submit(read_band_to_grid, s["assets"]["scl"], T, W, H, res, Resampling.nearest)
                bands = {b: f.result() for b, f in fb.items()}; scl = fs.result().astype(np.uint8)
            # quality mask expects the full band stack order; build a minimal one
            stack = np.zeros((len(C.BAND_ORDER), H, W), np.uint16)
            for b, arr in bands.items():
                stack[C.BAND_ORDER.index(b)] = arr
            for b in C.BAND_ORDER:            # fill unread bands so the "nodata" test only looks at what we read
                if b not in bands: stack[C.BAND_ORDER.index(b)] = np.where(bands["red"] > 0, 1, 0)
            valid = scene_quality_mask(stack, scl) | (np.isin(scl, [6]) & (bands["red"] > 0))   # keep water too
            take = valid & ~filled
            for i, b in enumerate(("red", "green", "blue")):
                rgb[i][take] = bands[b][take]
            filled |= take
            print(f"{s['id']}: filled {filled.mean():.3f}", flush=True)
            if filled.mean() > 0.995:
                break
        if filled.mean() > 0.995:
            break
    # stretch: reflectance 0..0.30 -> gamma 0.65, then unsharp mask for crispness
    img = np.clip(rgb.astype(np.float32) / C.BOA_SCALE / 0.30, 0, 1) ** 0.65
    for i in range(3):
        blur = ndimage.gaussian_filter(img[i], 1.2)
        img[i] = np.clip(img[i] + 0.6 * (img[i] - blur), 0, 1)
    img8 = (img * 255).astype(np.uint8); alpha = (filled * 255).astype(np.uint8)
    # reproject to Web Mercator at ~10 m ground resolution
    dst_crs = "EPSG:3857"
    dst_T, dw, dh = calculate_default_transform(C.UTM, dst_crs, W, H, left=T.c, bottom=T.f + H * T.e, right=T.c + W * T.a, top=T.f)
    merc = np.zeros((4, dh, dw), np.uint8)
    for i in range(3):
        reproject(img8[i], merc[i], src_transform=T, src_crs=C.UTM, dst_transform=dst_T, dst_crs=dst_crs, resampling=Resampling.bilinear)
    reproject(alpha, merc[3], src_transform=T, src_crs=C.UTM, dst_transform=dst_T, dst_crs=dst_crs, resampling=Resampling.nearest)
    from pyproj import Transformer
    to_ll = Transformer.from_crs(dst_crs, 4326, always_xy=True)
    tiles = []
    for r0 in range(0, dh, a.tile):
        for c0 in range(0, dw, a.tile):
            r1, c1 = min(r0 + a.tile, dh), min(c0 + a.tile, dw)
            t = merc[:, r0:r1, c0:c1]
            if t[3].max() == 0:
                continue
            x0, y1 = dst_T * (c0, r0); x1, y0 = dst_T * (c1, r1)
            w_, s_ = to_ll.transform(x0, y0); e_, n_ = to_ll.transform(x1, y1)
            name = f"tile_{r0 // a.tile}_{c0 // a.tile}.webp"
            Image.fromarray(np.moveaxis(t, 0, -1), "RGBA").save(os.path.join(out, name), "WEBP", quality=78, method=4)
            tiles.append(dict(file=name, bounds=[[s_, w_], [n_, e_]]))
    json.dump(dict(res_m=res, width=dw, height=dh, tiles=tiles, dates=dates), open(os.path.join(out, "base.json"), "w"))
    print(f"base: {dw}x{dh} px, {len(tiles)} tiles, {sum(os.path.getsize(os.path.join(out, t['file'])) for t in tiles) / 1e6:.1f} MB", flush=True)

if __name__ == "__main__":
    main()
