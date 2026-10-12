"""Drone-orthoimage reference for validating the colour layer.

Classifies every 10 cm pixel of an RGB ortho as shadow / non-vegetation / green leaf /
turned leaf (yellow-orange-red) / brown, then aggregates the fractions onto the 10 m
analysis grid of a region so they can be compared cell by cell with the satellite layers.

Usage: python ortho_reference.py ORTHO.tif REGION OUTDIR
Writes OUTDIR/ortho_cells.npz (per-cell fractions on the region grid), preview images,
and OUTDIR/ortho_meta.json (thresholds and footprint).
"""
import os, sys, json
import numpy as np, rasterio
from rasterio.windows import Window
from PIL import Image
import config as C
from ingest_s2 import region_grid

# Hue classes in degrees (HSV from camera RGB): leaf green sits ~70-150, yellow ~45-70, orange/red < 45 (and > 330)
H_TURN_MAX, H_GREEN_MIN, H_GREEN_MAX = 62, 62, 165
S_MIN, V_SHADOW, V_MAX_LEAF = 0.18, 0.16, 0.97

def classify(rgb):
    """rgb uint8 [H,W,3] -> class uint8: 0 nodata, 1 shadow, 2 non-veg, 3 green leaf, 4 turned leaf, 5 brown/bare."""
    f = rgb.astype(np.float32) / 255.0
    r, g, b = f[..., 0], f[..., 1], f[..., 2]
    mx = f.max(axis=-1); mn = f.min(axis=-1); d = mx - mn + 1e-6
    hue = np.where(mx == r, ((g - b) / d) % 6, np.where(mx == g, (b - r) / d + 2, (r - g) / d + 4)) * 60
    sat = d / (mx + 1e-6); val = mx
    exg = 2 * g - r - b                                       # excess green: vegetation vs built/soil
    out = np.full(r.shape, 2, np.uint8)
    nodata = (rgb.max(axis=-1) == 0); shadow = val < V_SHADOW
    green = (hue >= H_GREEN_MIN) & (hue < H_GREEN_MAX) & (sat >= S_MIN) & (exg > 0.04)
    turned = ((hue < H_TURN_MAX) | (hue > 330)) & (sat >= 0.30) & (val > 0.25) & (r > b * 1.25)
    brown = ((hue < 45) & (sat >= 0.18) & (sat < 0.30) & (val < 0.6)) | ((hue < 45) & (val < 0.35) & (sat >= 0.18))
    out[brown] = 5; out[green] = 3; out[turned] = 4
    out[shadow] = 1; out[nodata] = 0
    return out, hue

def main(ortho, region, outdir):
    os.makedirs(outdir, exist_ok=True)
    T, W, H = region_grid(region); res = C.REGIONS[region]["res"]
    counts = np.zeros((6, H, W), np.int64); hue_sum = np.zeros((H, W), np.float64); hue_n = np.zeros((H, W), np.int64)
    with rasterio.open(ortho) as ds:
        assert ds.crs.to_string() == C.UTM, ds.crs
        px = ds.res[0]; b = ds.bounds
        # cell indices of the ortho footprint on the region grid
        c0 = int((b.left - T.c) / res); c1 = int(np.ceil((b.right - T.c) / res)); r0 = int((T.f - b.top) / res); r1 = int(np.ceil((T.f - b.bottom) / res))
        step = 400   # ortho rows per block (10 cm -> 40 m)
        prev = np.zeros((ds.height // 10, ds.width // 10, 3), np.uint8); cls_prev = np.zeros((ds.height // 10, ds.width // 10), np.uint8)
        for y0 in range(0, ds.height, step):
            h = min(step, ds.height - y0)
            rgb = np.moveaxis(ds.read([1, 2, 3], window=Window(0, y0, ds.width, h)), 0, -1)
            cls, hue = classify(rgb)
            # preview at 1 m
            n = (y0 + h) // 10 - y0 // 10
            prev[y0 // 10:y0 // 10 + n] = rgb[::10, ::10][:n, :prev.shape[1]]
            cls_prev[y0 // 10:y0 // 10 + n] = cls[::10, ::10][:n, :prev.shape[1]]
            # aggregate to grid cells: map each ortho pixel to a cell
            ys = b.top - (y0 + np.arange(h) + 0.5) * px; xs = b.left + (np.arange(ds.width) + 0.5) * px
            rr = ((T.f - ys) / res).astype(int); cc = ((xs - T.c) / res).astype(int)
            ok_r = (rr >= 0) & (rr < H); ok_c = (cc >= 0) & (cc < W)
            R = np.repeat(rr[:, None], ds.width, axis=1); Cc = np.repeat(cc[None, :], h, axis=0); okm = ok_r[:, None] & ok_c[None, :]
            flat = (R * W + Cc)[okm]
            for k in range(6):
                sel = (cls == k)[okm]
                counts[k].reshape(-1)[:] += np.bincount(flat[sel], minlength=H * W)
            leaf = ((cls == 3) | (cls == 4))[okm]
            hue_sum.reshape(-1)[:] += np.bincount(flat[leaf], weights=hue[okm][leaf], minlength=H * W)
            hue_n.reshape(-1)[:] += np.bincount(flat[leaf], minlength=H * W)
    total = counts.sum(axis=0); valid = counts[1:].sum(axis=0)
    with np.errstate(all="ignore"):
        covered = total > 0.9 * (res / px) ** 2         # cell fully inside the ortho
        veg = counts[3] + counts[4]
        out = dict(covered=covered, veg_fraction=np.where(valid > 0, veg / np.maximum(valid, 1), np.nan),
                   turned_fraction=np.where(veg > 0, counts[4] / np.maximum(veg, 1), np.nan),
                   brown_fraction=np.where(valid > 0, counts[5] / np.maximum(valid, 1), np.nan),
                   shadow_fraction=np.where(total > 0, counts[1] / np.maximum(total, 1), np.nan),
                   leaf_hue=np.where(hue_n > 0, hue_sum / np.maximum(hue_n, 1), np.nan),
                   counts=counts, transform=np.array(list(T)[:6]), res=res)
    np.savez_compressed(os.path.join(outdir, "ortho_cells.npz"), **out)
    Image.fromarray(prev).save(os.path.join(outdir, "ortho_preview_1m.jpg"), quality=85)
    pal = np.array([[0, 0, 0], [40, 40, 60], [150, 150, 150], [26, 128, 38], [250, 140, 13], [120, 85, 60]], np.uint8)
    Image.fromarray(pal[cls_prev]).save(os.path.join(outdir, "ortho_classes_1m.png"))
    cov = covered.sum(); vf = out["veg_fraction"][covered]; tf = out["turned_fraction"][covered]
    meta = dict(ortho=os.path.basename(ortho), region=region, cells_covered=int(cov), cell_rows=[r0, r1], cell_cols=[c0, c1],
                thresholds=dict(H_TURN_MAX=H_TURN_MAX, H_GREEN_MIN=H_GREEN_MIN, H_GREEN_MAX=H_GREEN_MAX, S_MIN=S_MIN, V_SHADOW=V_SHADOW),
                class_share={n: round(float(counts[i][covered].sum() / max(total[covered].sum(), 1)), 3) for i, n in enumerate(["nodata", "shadow", "nonveg", "green", "turned", "brown"])},
                veg_fraction_median=round(float(np.nanmedian(vf)), 3), turned_fraction_pct=[round(float(np.nanpercentile(tf, q)), 3) for q in (10, 50, 90)])
    json.dump(meta, open(os.path.join(outdir, "ortho_meta.json"), "w"), indent=1)
    print(json.dumps(meta, indent=1))

if __name__ == "__main__":
    main(*sys.argv[1:4])
