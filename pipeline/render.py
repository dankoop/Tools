"""Render analysis rasters to Web-Mercator PNG/JPEG overlays for the Leaflet map.

Usage: python render.py REGION ANALYSIS.npz OUTDIR [--max-px 3000]
Writes OUTDIR/<region>/<layer>.png|jpg and OUTDIR/<region>/layers.json (bounds + legend info).
"""
import os, sys, json, argparse
import numpy as np, rasterio
from rasterio.warp import reproject, Resampling, calculate_default_transform
from rasterio.transform import Affine
from PIL import Image
import config as C

RAMP_CS = [(0.10, 0.50, 0.15), (0.95, 0.85, 0.10), (0.98, 0.55, 0.05), (0.85, 0.10, 0.05)]      # green→yellow→orange→red
STAGE_COL = {1: (26, 128, 38), 2: (242, 217, 26), 3: (250, 140, 13), 4: (120, 85, 60)}   # same hues as the colour ramp
# Non-colour diagnostics use cool / neutral hues so orange and yellow only ever mean fall colour.
RAMP_DIV = [(0.36, 0.23, 0.54), (0.90, 0.90, 0.90), (0.12, 0.37, 0.42)]                           # purple→grey→teal
RAMP_TEAL = [(0.86, 0.90, 0.91), (0.44, 0.71, 0.76), (0.12, 0.37, 0.42), (0.03, 0.18, 0.23)]        # light→dark teal
RAMP_GREY = [(0.86, 0.86, 0.86), (0.35, 0.35, 0.35)]
# days-to-peak: past = grey, now = darkest teal, later = progressively lighter
PEAK_STOPS = [(0.0, (0.72, 0.72, 0.72)), (0.30, (0.55, 0.55, 0.55)), (0.375, (0.03, 0.20, 0.26)), (0.50, (0.15, 0.45, 0.52)),
              (0.75, (0.44, 0.71, 0.76)), (1.0, (0.84, 0.91, 0.93))]

def ramp_stops(x, stops):
    x = np.clip(x, 0, 1); out = np.zeros(x.shape + (3,), np.float32)
    for (p0, c0), (p1, c1) in zip(stops[:-1], stops[1:]):
        m = (x >= p0) & (x <= p1); f = ((x - p0) / max(p1 - p0, 1e-6))[..., None]
        out[m] = (np.array(c0) * (1 - f) + np.array(c1) * f)[m]
    return out * 255

def ramp(x, ramp_cols):
    x = np.clip(x, 0, 1) * (len(ramp_cols) - 1); i = np.clip(x.astype(int), 0, len(ramp_cols) - 2); f = (x - i)[..., None]
    rc = np.array(ramp_cols)
    return (rc[i] * (1 - f) + rc[i + 1] * f) * 255

def warp(a, src_T, src_crs, dst_T, dst_crs, shape, resampling, nodata=0):
    out = np.full(shape, nodata, a.dtype)
    reproject(a, out, src_transform=src_T, src_crs=src_crs, dst_transform=dst_T, dst_crs=dst_crs,
              resampling=resampling, src_nodata=nodata, dst_nodata=nodata)
    return out

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("region"); ap.add_argument("npz"); ap.add_argument("outdir")
    ap.add_argument("--max-px", type=int, default=3200)
    a = ap.parse_args()
    z = np.load(a.npz); T = Affine(*z["transform"][:6]); res = int(z["res"])
    H, W = z["veg"].shape
    out = os.path.join(a.outdir, a.region); os.makedirs(out, exist_ok=True)
    dst_crs = "EPSG:3857"
    dst_T, dw, dh = calculate_default_transform(C.UTM, dst_crs, W, H, left=T.c, bottom=T.f + H * T.e, right=T.c + W * T.a, top=T.f)
    scale = max(dw, dh) / a.max_px
    if scale > 1:
        dst_T = Affine(dst_T.a * scale, 0, dst_T.c, 0, dst_T.e * scale, dst_T.f); dw, dh = int(dw / scale), int(dh / scale)
    shape = (dh, dw)
    bounds = rasterio.transform.array_bounds(dh, dw, dst_T)   # (west, south, east, north) in 3857
    from pyproj import Transformer
    to_ll = Transformer.from_crs(dst_crs, 4326, always_xy=True)
    w, s = to_ll.transform(bounds[0], bounds[1]); e, n = to_ll.transform(bounds[2], bounds[3])
    layers = {"bounds": [[s, w], [n, e]], "width": dw, "height": dh, "layers": {}}
    veg = warp(z["veg"].astype(np.uint8), T, C.UTM, dst_T, dst_crs, shape, Resampling.nearest).astype(bool)
    region_mask = warp((z["sea_mask"] if a.region == "seattle" else z["kc_mask"]).astype(np.uint8), T, C.UTM, dst_T, dst_crs, shape, Resampling.nearest).astype(bool)
    veg &= region_mask

    def save_rgba(name, rgb, alpha, fmt="png"):
        img = np.dstack([np.clip(rgb, 0, 255).astype(np.uint8), (alpha * 255).astype(np.uint8)])
        p = os.path.join(out, f"{name}.{fmt}")
        if fmt == "png":
            Image.fromarray(img, "RGBA").quantize(colors=255, method=Image.Quantize.FASTOCTREE).save(p, optimize=True)
        else:
            Image.fromarray(img, "RGBA").save(p, "WEBP", quality=80, method=4)
        return os.path.basename(p)

    def q(a8, lo, hi):   # undo q8 quantisation
        return np.where(a8 > 0, (a8.astype(np.float32) - 1) / 254 * (hi - lo) + lo, np.nan)

    def continuous(name, key, lo, hi, ramp_cols, title, labels, mask_veg=True, alpha_from=None):
        a8 = warp(z[key], T, C.UTM, dst_T, dst_crs, shape, Resampling.average)
        v = q(a8, lo, hi); m = np.isfinite(v) & (veg if mask_veg else True)
        x = (v - lo) / (hi - lo)
        rgb = ramp(np.nan_to_num(x), ramp_cols)
        alpha = m.astype(np.float32) * (0.92 if alpha_from is None else np.clip(alpha_from(x), 0.15, 0.95))
        layers["layers"][name] = dict(file=save_rgba(name, rgb, alpha), title=title, legend=dict(type="ramp", colors=[[int(c * 255) for c in rc] for rc in ramp_cols], labels=labels))

    continuous("cs_now", "cs_now", 0, 1, RAMP_CS, "Fall colour now", ["green", "turning", "peak", "deep red"], alpha_from=lambda x: 0.25 + 0.75 * x)
    continuous("dhue", "dhue", -10, 60, RAMP_CS, "Hue shift since early Sept (°)", ["0°", "+20°", "+40°", "+60°"])
    continuous("dndvi", "dndvi", -0.6, 0.2, RAMP_DIV, "Greenness change since early Sept (NDVI)", ["-0.6 (leaf-off)", "0", "+0.2"])
    continuous("velocity", "velocity", -0.05, 0.15, RAMP_TEAL, "Change velocity (senescence index per week)", ["steady", "slow", "fast", "very fast"])
    continuous("decid", "decid", 0, 1, RAMP_TEAL[:3], "Deciduous likelihood", ["evergreen", "mixed", "deciduous"])
    continuous("confidence", "confidence", 0, 1, RAMP_GREY, "Prediction confidence", ["low", "high"])
    # peak timing classes, chronological: past (leaf-off brown) | now (darkest teal) | next 2 wk | 2-4 wk | later (lightest)
    TIMING_COL = {1: (120, 85, 60), 2: (8, 51, 66), 3: (38, 115, 130), 4: (112, 181, 194), 5: (205, 226, 230)}
    tm = warp(z["timing"], T, C.UTM, dst_T, dst_crs, shape, Resampling.mode); tm[~region_mask] = 0
    rgb = np.zeros(shape + (3,), np.float32)
    for k_, col in TIMING_COL.items():
        rgb[tm == k_] = col
    layers["layers"]["timing"] = dict(file=save_rgba("timing", rgb, (tm > 0) * 0.88), title="Predicted peak timing",
        legend=dict(type="classes", items=[["past peak", TIMING_COL[1]], ["peak now", TIMING_COL[2]], ["next 2 weeks", TIMING_COL[3]], ["2-4 weeks", TIMING_COL[4]], ["later", TIMING_COL[5]]]))
    conf = q(warp(z["confidence"], T, C.UTM, dst_T, dst_crs, shape, Resampling.average), 0, 1); conf = np.nan_to_num(conf, nan=0.3)
    confw = 0.45 + 0.55 * np.clip(conf, 0, 1)                   # forecast layers fade where the prediction is weak
    # stage (now)
    def stage_layer(name, st8):
        st = warp(st8, T, C.UTM, dst_T, dst_crs, shape, Resampling.mode); st[~region_mask] = 0
        rgb = np.zeros(shape + (3,), np.float32)
        for k_, col in STAGE_COL.items():
            rgb[st == k_] = col
        return save_rgba(name, rgb, (st > 0) * 0.85)
    layers["layers"]["stage"] = dict(file=stage_layer("stage", z["stage"]), title="Colour stage now",
        legend=dict(type="classes", items=[["green / evergreen", STAGE_COL[1]], ["turning", STAGE_COL[2]], ["peak colour", STAGE_COL[3]], ["past peak / leaf-off", STAGE_COL[4]]]))
    # weekly layers
    weeks = []
    for lab in z["bins"]:
        lab = str(lab); slug = lab.replace(" ", "_").replace("-", "_")
        cs8 = warp(z[f"wk|{lab}|cs"], T, C.UTM, dst_T, dst_crs, shape, Resampling.average); v = q(cs8, 0, 1); m = np.isfinite(v) & veg
        f_cs = save_rgba(f"wk_{slug}_cs", ramp(np.nan_to_num(v), RAMP_CS), m * np.clip(0.25 + 0.75 * np.nan_to_num(v), 0, 0.95))
        f_st = stage_layer(f"wk_{slug}_stage", z[f"wk|{lab}|stage"])
        rgb = np.stack([warp(z[f"wk|{lab}|rgb"][i], T, C.UTM, dst_T, dst_crs, shape, Resampling.average) for i in range(3)], axis=-1)
        f_rgb = save_rgba(f"wk_{slug}_rgb", rgb, ((rgb.max(axis=-1) > 0) & region_mask) * 1.0, fmt="webp")
        cov = float((z[f"wk|{lab}|nobs"] > 0)[z["veg"]].mean()) if z["veg"].any() else 0   # share of vegetation with a clear pass that week, on the native grid
        weeks.append(dict(label=lab, cs=f_cs, stage=f_st, rgb=f_rgb, coverage=round(cov, 3)))
    layers["weeks"] = weeks
    fcs = []
    for rec in z["fc_bins"]:
        lab, a0, a1 = str(rec).split("|"); slug = lab.replace(" ", "_").replace("-", "_")
        cs8 = warp(z[f"fc|{lab}|cs"], T, C.UTM, dst_T, dst_crs, shape, Resampling.average); v = q(cs8, 0, 1); m = np.isfinite(v) & veg
        f_cs = save_rgba(f"fc_{slug}_cs", ramp(np.nan_to_num(v), RAMP_CS), m * np.clip(0.25 + 0.75 * np.nan_to_num(v), 0, 0.95) * confw)
        st = warp(z[f"fc|{lab}|stage"], T, C.UTM, dst_T, dst_crs, shape, Resampling.mode); st[~region_mask] = 0
        rgb = np.zeros(shape + (3,), np.float32)
        for k_, col in STAGE_COL.items():
            rgb[st == k_] = col
        f_st = save_rgba(f"fc_{slug}_stage", rgb, (st > 0) * 0.85 * confw)
        fcs.append(dict(label=lab, start=a0, end=a1, cs=f_cs, stage=f_st))
    layers["forecast"] = fcs
    json.dump(layers, open(os.path.join(out, "layers.json"), "w"), indent=1)
    print(a.region, "rendered", dw, "x", dh, "->", out, [f"{k}: {os.path.getsize(os.path.join(out, v['file'])) // 1024} KB" for k, v in layers["layers"].items()][:4])

if __name__ == "__main__":
    main()
