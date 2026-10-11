"""Fall-colour analysis on an ingested Sentinel-2 stack.

Steps (all per region grid):
  1. Per-scene spectral indices from the 10 L2A bands.
  2. Weekly cloud-free composites (median) + an early-September baseline composite.
  3. "Now" layer = latest clear observation per pixel (with its date) blended with the
     latest weekly composite where single-date coverage is missing.
  4. Colour score, stage classification, change vs baseline, velocity.
  5. Deciduous likelihood from the 2023-2025 baseline years (amplitude of autumn change).
  6. Phenology: climatological peak timing per pixel from baseline years, phase offset
     of this year, forecast-adjusted days-to-peak.
  7. Hotspot clustering and destination scoring -> GeoJSON.

Usage:
  python analyze.py composite  STACKDIR OUT.npz [--bins current|baseline]
  python analyze.py phenology  OUT.npz COMPOSITE_2023.npz COMPOSITE_2024.npz ...
  python analyze.py region     REGION STACKDIR COMPOSITE.npz PHENO.npz OUTDIR
"""
import os, sys, json, glob, argparse, datetime as dt, math
import numpy as np
from scipy import ndimage
import config as C

CACHE = C.SCRATCH
TODAY = dt.date.today()
BO = {b: i for i, b in enumerate(C.BAND_ORDER)}

# ----------------------------------------------------------------------------- indices
def reflectance(bands):
    return (bands.astype(np.float32) - C.BOA_OFFSET) / C.BOA_SCALE

def hue_chroma(r, g, b):
    mx = np.maximum(np.maximum(r, g), b); mn = np.minimum(np.minimum(r, g), b); d = mx - mn + 1e-6
    h = np.where(mx == r, ((g - b) / d) % 6, np.where(mx == g, (b - r) / d + 2, (r - g) / d + 4)) * 60
    return h.astype(np.float32), ((mx - mn) / (mx + 1e-6)).astype(np.float32)

INDEX_NAMES = ["ndvi", "psri", "cire", "ndmi", "hue", "chroma", "r", "g", "b"]

def indices(bands, valid):
    """Return dict of float32 arrays (NaN where invalid)."""
    p = reflectance(bands)
    blue, green, red = p[BO["blue"]], p[BO["green"]], p[BO["red"]]
    re1, re2, nir, nir8 = p[BO["rededge1"]], p[BO["rededge2"]], p[BO["nir"]], p[BO["nir08"]]
    swir = p[BO["swir16"]]
    out = {}
    with np.errstate(all="ignore"):
        out["ndvi"] = (nir - red) / (nir + red + 1e-6)
        out["psri"] = (red - blue) / (re2 + 1e-6)                 # plant senescence reflectance index
        out["cire"] = nir8 / (re1 + 1e-6) - 1                     # red-edge chlorophyll index
        out["ndmi"] = (nir8 - swir) / (nir8 + swir + 1e-6)        # leaf moisture / leaf-off
        out["hue"], out["chroma"] = hue_chroma(red, green, blue)
        out["r"], out["g"], out["b"] = red, green, blue
    # thin cloud / haze over vegetation: green leaves are dark in blue (< ~0.06); haze lifts it
    with np.errstate(all="ignore"):
        haze = (blue > 0.09) & (out["ndvi"] > 0.25)
        haze = ndimage.binary_dilation(haze, iterations=1)
    valid = valid & ~haze
    for k in out:
        out[k] = np.where(valid, out[k], np.nan).astype(np.float32)
    return out

# ----------------------------------------------------------------------------- composites
def current_bins(dates):
    """Weekly bins for the current season: [Sep1-14 baseline], then 7-day bins to the last date."""
    bins = [("Sep 1-14", "09-01", "09-14")]
    start = dt.date(2000, 9, 15)
    last = max(dates)
    while True:
        end = start + dt.timedelta(days=6)
        lab = f"{start.strftime('%b')} {start.day}-{end.day}" if start.month == end.month else f"{start.strftime('%b %-d')}-{end.strftime('%b %-d')}"
        bins.append((lab, start.strftime("%m-%d"), end.strftime("%m-%d")))
        if end.strftime("%m-%d") >= last[5:]:
            break
        start = end + dt.timedelta(days=1)
    return bins

def baseline_bins():
    bins = [("Sep 1-14", "09-01", "09-14")]
    start = dt.date(2000, 9, 15)
    while start <= dt.date(2000, 11, 30):
        end = start + dt.timedelta(days=6)
        bins.append((f"{start.strftime('%m-%d')}", start.strftime("%m-%d"), end.strftime("%m-%d")))
        start = end + dt.timedelta(days=1)
    return bins

def composite(stackdir, out_path, mode):
    files = sorted(glob.glob(os.path.join(stackdir, "*.npz")))
    dates = [os.path.basename(f)[:10] for f in files]
    bins = current_bins(dates) if mode == "current" else baseline_bins()
    grid = json.load(open(os.path.join(stackdir, "grid.json")))
    H, W = grid["height"], grid["width"]
    res = {}
    # latest-valid "now" layer and per-date accumulators for velocity (last 16 days)
    latest = {k: np.full((H, W), np.nan, np.float32) for k in INDEX_NAMES}
    latest_doy = np.zeros((H, W), np.int16)
    nobs_total = np.zeros((H, W), np.uint8)
    last_date = dt.date.fromisoformat(dates[-1])
    vel = {k: np.zeros((H, W), np.float64) for k in ("n", "st", "sy", "sty", "stt")}
    for lab, b0, b1 in bins:
        sel = [f for f, d in zip(files, dates) if b0 <= d[5:] <= b1]
        if not sel:
            res[lab] = None; continue
        acc = {k: [] for k in INDEX_NAMES}
        for f in sel:
            z = np.load(f); v = z["valid"]
            if v.sum() == 0:
                continue
            ix = indices(z["bands"], v)
            for k in INDEX_NAMES:
                acc[k].append(ix[k])
            if mode == "current":
                d = dt.date.fromisoformat(os.path.basename(f)[:10])
                for k in INDEX_NAMES:
                    latest[k][v] = ix[k][v]
                latest_doy[v] = d.timetuple().tm_yday
                nobs_total += v
                t = (d - last_date).days
                if t >= -16:
                    y = ix["psri"]; m = v & np.isfinite(y)
                    vel["n"][m] += 1; vel["st"][m] += t; vel["sy"][m] += y[m]; vel["sty"][m] += t * y[m]; vel["stt"][m] += t * t
        if not acc["ndvi"]:
            res[lab] = None; continue
        with np.errstate(all="ignore"):
            comp = {k: np.nanmedian(np.stack(acc[k]), axis=0).astype(np.float32) for k in INDEX_NAMES}
        comp["nobs"] = np.sum([np.isfinite(a) for a in acc["ndvi"]], axis=0).astype(np.uint8)
        res[lab] = comp
        print(f"  bin {lab}: {len(sel)} scenes, coverage {np.isfinite(comp['ndvi']).mean():.3f}", flush=True)
    out = {"bins": np.array([b[0] for b in bins]), "bin_ranges": np.array([f"{b[1]}|{b[2]}" for b in bins]),
           "dates": np.array(dates), "transform": np.array(grid["transform"]), "res": grid["res"]}
    for lab, comp in res.items():
        if comp is None:
            continue
        for k, a in comp.items():
            out[f"{lab}|{k}"] = a
    if mode == "current":
        for k in INDEX_NAMES:
            out[f"latest|{k}"] = latest[k]
        out["latest|doy"] = latest_doy; out["nobs_total"] = nobs_total
        with np.errstate(all="ignore"):
            n = vel["n"]; slope = (n * vel["sty"] - vel["st"] * vel["sy"]) / (n * vel["stt"] - vel["st"] ** 2)
        out["velocity_psri"] = np.where(n >= 2, slope, np.nan).astype(np.float32)
    np.savez_compressed(out_path, **out)
    print("wrote", out_path)

# ----------------------------------------------------------------------------- colour score
def colour_score(hue, chroma, psri, ndvi, hue0, psri0):
    """0..1 'how much fall colour' relative to a green baseline (hue0/psri0)."""
    with np.errstate(all="ignore"):
        hs = np.clip((hue0 - hue) / 50.0, 0, 1)            # hue shift green -> yellow/orange/red
        ps = np.clip((psri - psri0) / 0.15, 0, 1)          # senescence index rise
        leaves = np.clip((ndvi - 0.15) / 0.15, 0, 1)       # leaves still present
        cs = (0.6 * hs + 0.4 * ps) * leaves
    return cs.astype(np.float32)

def load_bin(z, lab, k):
    key = f"{lab}|{k}"
    return z[key] if key in z.files else None

# ----------------------------------------------------------------------------- phenology
def phenology(out_path, comp_paths):
    """Climatology from baseline years: per pixel colour-score curve by weekly bin, peak DOY,
    amplitude, autumn NDVI drop. All on the baseline (80 m) grid."""
    curves, peaks, amps, drops, leafoff = [], [], [], [], []
    bins = None
    for p in comp_paths:
        z = np.load(p)
        bins = list(z["bins"]); base = bins[0]
        hue0, psri0, ndvi0 = z[f"{base}|hue"], z[f"{base}|psri"], z[f"{base}|ndvi"]
        cs_w, ndvi_w = [], []
        for lab in bins[1:]:
            h = load_bin(z, lab, "hue")
            if h is None:
                cs_w.append(np.full(hue0.shape, np.nan, np.float32)); ndvi_w.append(np.full(hue0.shape, np.nan, np.float32)); continue
            cs_w.append(colour_score(h, z[f"{lab}|chroma"], z[f"{lab}|psri"], z[f"{lab}|ndvi"], hue0, psri0))
            ndvi_w.append(z[f"{lab}|ndvi"])
        cs_w = np.stack(cs_w); ndvi_w = np.stack(ndvi_w)
        # fill gaps along time by linear interpolation
        cs_w = interp_time(cs_w); ndvi_w = interp_time(ndvi_w)
        amp = np.nanmax(cs_w, axis=0)
        pk = np.nanargmax(np.nan_to_num(cs_w, nan=-1), axis=0)
        mid = [dt.date(2001, int(b[:2]), int(b[3:5])) + dt.timedelta(days=3) for b in [r.split("|")[0] for r in z["bin_ranges"][1:]]]
        doys = np.array([m.timetuple().tm_yday for m in mid])
        peaks.append(doys[pk].astype(np.float32)); amps.append(amp)
        drop = ndvi0 - np.nanmin(ndvi_w, axis=0); drops.append(drop)
        curves.append(cs_w / np.maximum(amp, 0.05)[None])
        print(f"  {os.path.basename(p)}: median amp {np.nanmedian(amp):.3f}, median peak DOY {np.nanmedian(doys[pk]):.0f}")
    curve = np.nanmean(np.stack(curves), axis=0).astype(np.float32)
    np.savez_compressed(out_path, curve=curve, doys=doys, peak_doy=np.nanmedian(np.stack(peaks), axis=0).astype(np.float32),
                        amp=np.nanmedian(np.stack(amps), axis=0).astype(np.float32),
                        ndvi_drop=np.nanmedian(np.stack(drops), axis=0).astype(np.float32),
                        transform=np.load(comp_paths[0])["transform"], res=np.load(comp_paths[0])["res"])
    print("wrote", out_path)

def interp_time(a):
    """Linear interpolation of NaNs along axis 0 (small number of time steps)."""
    T = a.shape[0]; out = a.copy()
    idx = np.arange(T)
    valid = np.isfinite(a)
    # forward/backward fill using nearest valid index on each side
    prev = np.where(valid, idx[:, None, None], -1); prev = np.maximum.accumulate(prev, axis=0)
    nxt = np.where(valid, idx[:, None, None], T); nxt = np.minimum.accumulate(nxt[::-1], axis=0)[::-1]
    flat = a.reshape(T, -1); pf = prev.reshape(T, -1); nf = nxt.reshape(T, -1); of = out.reshape(T, -1)
    for t in range(T):
        m = ~np.isfinite(flat[t]) & (pf[t] >= 0) & (nf[t] < T)
        cols = np.where(m)[0]
        if cols.size:
            p = pf[t, cols]; n = nf[t, cols]
            w = (t - p) / np.maximum(n - p, 1)
            of[t, cols] = flat[p, cols] * (1 - w) + flat[n, cols] * w
        m2 = ~np.isfinite(flat[t]) & (pf[t] >= 0) & (nf[t] >= T); cols = np.where(m2)[0]
        if cols.size: of[t, cols] = flat[pf[t, cols], cols]
        m3 = ~np.isfinite(flat[t]) & (pf[t] < 0) & (nf[t] < T); cols = np.where(m3)[0]
        if cols.size: of[t, cols] = flat[nf[t, cols], cols]
    return out

def resample_to_grid(src, src_transform, src_res, dst_transform, dst_shape):
    """Nearest-neighbour resample of a (…, H, W) array between two aligned UTM grids."""
    H, W = dst_shape
    xs = dst_transform[2] + (np.arange(W) + 0.5) * dst_transform[0]
    ys = dst_transform[5] + (np.arange(H) + 0.5) * dst_transform[4]
    cols = np.clip(((xs - src_transform[2]) / src_transform[0]).astype(int), 0, src.shape[-1] - 1)
    rows = np.clip(((ys - src_transform[5]) / src_transform[4]).astype(int), 0, src.shape[-2] - 1)
    return src[..., rows[:, None], cols[None, :]]

# ----------------------------------------------------------------------------- region analysis
STAGE_NAMES = ["masked", "green", "turning", "peak", "past peak"]

def region(name, stackdir, comp_path, pheno_path, outdir):
    os.makedirs(outdir, exist_ok=True)
    z = np.load(comp_path); anc = np.load(os.path.join(CACHE, f"anc_{name}.npz"))
    ph = np.load(pheno_path); weather = json.load(open(os.path.join(CACHE, "weather", "weather.json")))
    bins = list(z["bins"]); base = bins[0]; T = z["transform"]; res = int(z["res"])
    H, W = z[f"{base}|ndvi"].shape
    lc, dem, kc, sea, parks = anc["landcover"], anc["dem"], anc["kc_mask"], anc["seattle_mask"], anc["park_mask"]
    hue0, psri0, ndvi0, chroma0 = z[f"{base}|hue"], z[f"{base}|psri"], z[f"{base}|ndvi"], z[f"{base}|chroma"]
    region_mask = kc if name == "county" else sea
    veg = (ndvi0 > 0.4) & ~np.isin(lc, [40, 70, 80]) & region_mask & np.isfinite(hue0)

    # -- phenology layers resampled to this grid
    pt = ph["transform"]
    amp = resample_to_grid(ph["amp"], pt, int(ph["res"]), T, (H, W))
    drop = resample_to_grid(ph["ndvi_drop"], pt, int(ph["res"]), T, (H, W))
    peak_doy = resample_to_grid(ph["peak_doy"], pt, int(ph["res"]), T, (H, W))
    curve = resample_to_grid(ph["curve"], pt, int(ph["res"]), T, (H, W))      # [Tb, H, W]
    doys = ph["doys"]

    # -- weekly colour score + stage
    weekly = {}
    for lab in bins[1:]:
        h = load_bin(z, lab, "hue")
        if h is None:
            continue
        cs = colour_score(h, z[f"{lab}|chroma"], z[f"{lab}|psri"], z[f"{lab}|ndvi"], hue0, psri0)
        weekly[lab] = dict(cs=cs, ndvi=z[f"{lab}|ndvi"], hue=h, nobs=z[f"{lab}|nobs"], r=z[f"{lab}|r"], g=z[f"{lab}|g"], b=z[f"{lab}|b"])
    # -- "now": latest observation, fall back to last two weekly composites
    now = {k: z[f"latest|{k}"].copy() for k in ("hue", "chroma", "psri", "ndvi", "ndmi", "cire")}
    now_doy = z["latest|doy"].astype(np.int16)
    stale = now_doy < (TODAY - dt.timedelta(days=12)).timetuple().tm_yday
    for lab in reversed(list(weekly)[-2:]):
        wk = weekly[lab]
        fill = (~np.isfinite(now["hue"]) | stale) & np.isfinite(wk["hue"])
        for k in ("hue", "ndvi"):
            now[k][fill] = wk[k][fill]
        now["psri"][fill] = z[f"{lab}|psri"][fill]; now["chroma"][fill] = z[f"{lab}|chroma"][fill]
        now["ndmi"][fill] = z[f"{lab}|ndmi"][fill]; now["cire"][fill] = z[f"{lab}|cire"][fill]
        stale[fill] = False
    cs_now = colour_score(now["hue"], now["chroma"], now["psri"], now["ndvi"], hue0, psri0)
    with np.errstate(all="ignore"):
        dndvi = now["ndvi"] - ndvi0; dhue = hue0 - now["hue"]; dpsri = now["psri"] - psri0
    # deciduous likelihood: prior-year amplitude / NDVI drop, boosted by this year's change
    decid = np.maximum(np.clip((amp - 0.12) / 0.35, 0, 1), np.clip((drop - 0.12) / 0.3, 0, 1))
    decid = np.where(cs_now > 0.35, np.maximum(decid, 0.8), decid)
    decid = np.where(veg, decid, 0).astype(np.float32)
    # stage
    stage = np.zeros((H, W), np.uint8)
    ok = veg & np.isfinite(cs_now)
    stage[ok] = 1
    stage[ok & (cs_now >= 0.15)] = 2
    stage[ok & (cs_now >= 0.45) & (now["ndvi"] >= 0.3)] = 3
    stage[ok & (dndvi <= -0.3) & (now["ndvi"] < 0.35)] = 4
    stage[veg & (decid < 0.25) & (cs_now < 0.15)] = 1

    # -- phenology: phase offset and days to peak
    this_doy = TODAY.timetuple().tm_yday
    with np.errstate(all="ignore"):
        p = np.clip(cs_now / np.maximum(amp, 0.15), 0, 1.0)
    # equivalent DOY = first bin where climatological normalised curve reaches p (interpolated)
    cur = np.nan_to_num(curve, nan=0.0)
    Tb = cur.shape[0]
    eq = np.full((H, W), np.nan, np.float32)
    reached = np.zeros((H, W), bool)
    for i in range(Tb - 1):
        c0, c1 = cur[i], cur[i + 1]
        hit = ~reached & (c1 >= p) & (c1 > c0)
        w = np.clip((p - c0) / np.maximum(c1 - c0, 1e-3), 0, 1)
        eq[hit] = (doys[i] + w * (doys[i + 1] - doys[i]))[hit]
        reached |= hit
    eq[~reached & (p < 0.05)] = doys[0] - 7           # nothing happened yet: before season start
    obs_doy = np.where(now_doy > 0, now_doy, this_doy).astype(np.float32)
    offset = eq - obs_doy                              # +: ahead of a normal year
    offset = np.clip(np.nan_to_num(offset, nan=0.0), -30, 30)
    pred_peak = peak_doy - offset
    days_to_peak = pred_peak - this_doy
    # weather adjustment: forecast cold nights accelerate; lapse-rate adjust by elevation
    gfs = weather["gfs"]["Seattle"]
    cold = sum(1 for d in gfs[:7] if d["tmin"] is not None and d["tmin"] <= 5)
    elev_adj = np.clip((dem - 50) / 1000.0 * 6.5, 0, 12)     # colder with elevation (deg C)
    cold_px = cold + np.where(elev_adj > 3, 2, np.where(elev_adj > 1.5, 1, 0))
    days_to_peak = days_to_peak - np.clip(0.5 * cold_px, 0, 4)
    days_to_peak = np.where(veg & (decid > 0.25), days_to_peak, np.nan).astype(np.float32)
    k_sm = 7 if res == 10 else 5
    num = ndimage.uniform_filter(np.nan_to_num(days_to_peak), k_sm); den = ndimage.uniform_filter(np.isfinite(days_to_peak).astype(np.float32), k_sm)
    with np.errstate(all="ignore"):
        days_to_peak = np.where(np.isfinite(days_to_peak), num / np.maximum(den, 1e-3), np.nan).astype(np.float32)
    days_to_peak = np.where(stage == 4, np.minimum(days_to_peak, -10), days_to_peak)
    conf = np.clip(1 - np.abs(offset) / 30, 0.2, 1) * np.clip(amp / 0.4, 0.3, 1)
    gust = max((d["gust_max_ms"] or 0) for d in gfs[:7]); rain = max(d["precip_mm"] for d in gfs[:7])
    leafdrop_risk = "high" if gust >= 18 or rain >= 25 else "moderate" if gust >= 12 or rain >= 12 else "low"

    # -- save rasters
    def q8(a, lo, hi):
        return np.where(np.isfinite(a), np.clip((a - lo) / (hi - lo) * 254 + 1, 1, 255), 0).astype(np.uint8)
    out = dict(transform=np.array(T), res=res, bins=np.array(list(weekly)), veg=veg, stage=stage,
               cs_now=q8(cs_now, 0, 1), dhue=q8(dhue, -10, 60), dndvi=q8(dndvi, -0.6, 0.2), dpsri=q8(dpsri, -0.1, 0.4),
               velocity=q8(z["velocity_psri"] * 7, -0.05, 0.15), days_to_peak=q8(days_to_peak, -21, 35),
               decid=q8(decid, 0, 1), obs_doy=now_doy, confidence=q8(conf, 0, 1), dem=dem.astype(np.float32),
               nobs_total=z["nobs_total"], latest_doy_max=int(now_doy.max()))
    for lab, wk in weekly.items():
        out[f"wk|{lab}|cs"] = q8(wk["cs"], 0, 1); out[f"wk|{lab}|nobs"] = wk["nobs"]
        st = np.zeros((H, W), np.uint8); okw = veg & np.isfinite(wk["cs"])
        st[okw] = 1; st[okw & (wk["cs"] >= 0.15)] = 2; st[okw & (wk["cs"] >= 0.45) & (wk["ndvi"] >= 0.3)] = 3
        st[okw & ((wk["ndvi"] - ndvi0) <= -0.3) & (wk["ndvi"] < 0.35)] = 4
        out[f"wk|{lab}|stage"] = st
        rgb = np.stack([wk["r"], wk["g"], wk["b"]])
        with np.errstate(all="ignore"):
            out[f"wk|{lab}|rgb"] = np.where(np.isfinite(rgb), np.clip((rgb / 0.25) ** 0.7 * 255, 1, 255), 0).astype(np.uint8)
    np.savez_compressed(os.path.join(outdir, f"analysis_{name}.npz"), **out)

    # -- stats
    vv = veg & (decid > 0.25)
    summary = dict(region=name, today=str(TODAY), latest_obs=str(dt.date(TODAY.year, 1, 1) + dt.timedelta(int(now_doy.max()) - 1)),
                   coverage_latest_week=float(np.isfinite(weekly[list(weekly)[-1]]["cs"])[veg].mean()),
                   stage_fraction={STAGE_NAMES[i]: float((stage[vv] == i).mean()) for i in range(1, 5)},
                   median_days_to_peak=float(np.nanmedian(days_to_peak[vv])),
                   forecast_cold_nights=cold, forecast_max_gust_ms=gust, forecast_max_rain_mm=rain, leafdrop_risk=leafdrop_risk,
                   weekly_mean_cs={lab: float(np.nanmean(wk["cs"][vv])) for lab, wk in weekly.items()},
                   weekly_stage_fraction={lab: {STAGE_NAMES[i]: float((out[f'wk|{lab}|stage'][vv] == i).mean()) for i in range(1, 5)} for lab in weekly})
    json.dump(summary, open(os.path.join(outdir, f"summary_{name}.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in summary.items() if k not in ("weekly_stage_fraction",)}, indent=1))

    clusters(name, outdir, T, res, veg, decid, cs_now, stage, days_to_peak, z["velocity_psri"], dem, parks, sea, lc, weekly, anc, summary)

# ----------------------------------------------------------------------------- clusters
def clusters(name, outdir, T, res, veg, decid, cs_now, stage, dtp, vel, dem, parks, sea, lc, weekly, anc, summary):
    from rasterio import features
    from rasterio.transform import Affine
    from shapely.geometry import shape, mapping, Point
    from shapely.ops import unary_union
    from pyproj import Transformer
    cell_m = 100 if name == "seattle" else 200
    k = cell_m // res
    H, W = veg.shape; Hc, Wc = H // k, W // k
    def block(a, fn):
        a = a[:Hc * k, :Wc * k].reshape(Hc, k, Wc, k)
        return fn(a)
    vv = (veg & (decid > 0.25)).astype(np.float32)
    nveg = block(vv, lambda a: a.sum(axis=(1, 3)))
    csw = np.where(vv > 0, np.nan_to_num(cs_now), 0)
    cs_cell = block(csw, lambda a: a.sum(axis=(1, 3))) / np.maximum(nveg, 1)
    frac = nveg / (k * k)
    score = cs_cell * np.clip(frac / 0.25, 0, 1)
    peak_cell = block(np.where(vv > 0, (stage == 3).astype(np.float32), 0), lambda a: a.sum(axis=(1, 3))) / np.maximum(nveg, 1)
    dtp_cell = block(np.where(vv > 0, np.nan_to_num(dtp, nan=0), 0), lambda a: a.sum(axis=(1, 3))) / np.maximum(nveg, 1)
    vel_cell = block(np.where(vv > 0, np.nan_to_num(vel, nan=0), 0), lambda a: a.sum(axis=(1, 3))) / np.maximum(nveg, 1)
    dem_cell = block(np.nan_to_num(dem), lambda a: a.mean(axis=(1, 3)))
    park_cell = block(parks.astype(np.float32), lambda a: a.mean(axis=(1, 3)))
    sea_cell = block(sea.astype(np.float32), lambda a: a.mean(axis=(1, 3)))
    built = block((lc == 50).astype(np.float32), lambda a: a.mean(axis=(1, 3)))
    near_built = ndimage.maximum_filter(built, size=int(1500 / cell_m) * 2 + 1) > 0.15
    thr = max(0.22, float(np.nanpercentile(score[frac > 0.2], 90)))
    hot = (score >= thr) & (frac >= 0.2)
    hot = ndimage.binary_opening(hot, iterations=1) | (score >= thr * 1.3) & (frac >= 0.2)
    lab, n = ndimage.label(hot, structure=np.ones((3, 3)))
    print(f"[{name}] hotspot threshold {thr:.3f}; {n} raw clusters")
    cell_T = Affine(T[0] * k, 0, T[2], 0, T[4] * k, T[5])
    to_ll = Transformer.from_crs(C.UTM, 4326, always_xy=True)
    import shapely
    places = json.load(open(os.path.join(os.path.dirname(__file__), "data", "places.json")))
    nh_names = list(anc["nhood_names"]); nh_id = anc["nhood_id"]
    nh_cell = block(nh_id, lambda a: a[:, 0, :, 0])
    place_pts = [(p["name"], Point(*Transformer.from_crs(4326, C.UTM, always_xy=True).transform(p["lon"], p["lat"]))) for p in places["places"] + places["destinations"]]
    feats = []
    min_cells = 4 if name == "seattle" else 3
    for i in range(1, n + 1):
        m = lab == i
        nc = int(m.sum())
        if nc < min_cells:
            continue
        geoms = [shape(g) for g, v in features.shapes(m.astype(np.uint8), mask=m, transform=cell_T)]
        poly = unary_union(geoms).simplify(cell_m * 0.6)
        cen = poly.centroid
        area_ha = nc * cell_m * cell_m / 1e4
        intensity = float(np.average(cs_cell[m], weights=nveg[m]))
        peakf = float(np.average(peak_cell[m], weights=nveg[m])); d2p = float(np.average(dtp_cell[m], weights=nveg[m]))
        velo = float(np.average(vel_cell[m], weights=nveg[m])) * 7
        elev = float(dem_cell[m].mean()); parkf = float(park_cell[m].mean()); seaf = float(sea_cell[m].mean())
        access = max(parkf, 0.6 * float(near_built[m].mean()))
        ids = nh_cell[m]; ids = ids[ids > 0]
        if ids.size and seaf > 0.3:
            nh = str(nh_names[int(np.bincount(ids).argmax()) - 1]); label = nh
        else:
            dists = [(cen.distance(pt), nm) for nm, pt in place_pts]; dmin, nm = min(dists); label = f"near {nm}" if dmin > 1500 else nm
        timing = 1.0 if abs(d2p) <= 7 else max(0.0, 1 - (abs(d2p) - 7) / 14) if d2p > 0 else max(0.1, 1 - (abs(d2p) - 7) / 10)
        visit = 100 * (0.45 * min(intensity / 0.6, 1) + 0.2 * timing + 0.2 * access + 0.15 * min(math.log1p(area_ha) / math.log1p(400), 1))
        series = {wl: float(np.nanmean(np.where(vv[:Hc * k, :Wc * k] > 0, weekly[wl]["cs"][:Hc * k, :Wc * k], np.nan)[np.kron(m, np.ones((k, k), bool))])) for wl in weekly}
        poly_ll = shapely.transform(poly, lambda xy: np.column_stack(to_ll.transform(xy[:, 0], xy[:, 1])))
        feats.append(dict(type="Feature", geometry=mapping(poly_ll), properties=dict(
            label=label, area_ha=round(area_ha, 1), intensity=round(intensity, 3), peak_fraction=round(peakf, 3),
            days_to_peak=round(d2p, 1), velocity_wk=round(velo, 3), elev_m=round(elev), park_fraction=round(parkf, 2),
            in_seattle=seaf > 0.5, access=round(access, 2), visit_score=round(visit, 1), region=name,
            series={k_: (None if np.isnan(v) else round(v, 3)) for k_, v in series.items()},
            lon=round(float(poly_ll.centroid.x), 5), lat=round(float(poly_ll.centroid.y), 5))))
    feats.sort(key=lambda f: -f["properties"]["visit_score"])
    for r, f in enumerate(feats, 1):
        f["properties"]["rank"] = r
    json.dump(dict(type="FeatureCollection", features=feats), open(os.path.join(outdir, f"clusters_{name}.geojson"), "w"))
    print(f"[{name}] {len(feats)} clusters; top:", [(f['properties']['label'], f['properties']['visit_score']) for f in feats[:8]])

    # destinations
    if name == "county" or True:
        tr = Transformer.from_crs(4326, C.UTM, always_xy=True)
        dest = []
        for d in places["destinations"]:
            x, y = tr.transform(d["lon"], d["lat"]); c = int((x - T[2]) / T[0]); r = int((y - T[5]) / T[4]); rad = int(d["radius_m"] / res)
            if not (0 <= r < H and 0 <= c < W):
                continue
            if name == "seattle" and not sea[r, c]:
                continue
            if name == "county" and sea[r, c]:
                continue   # scored on the 10 m Seattle grid instead
            r0, r1, c0, c1 = max(r - rad, 0), min(r + rad, H), max(c - rad, 0), min(c + rad, W)
            yy, xx = np.ogrid[r0:r1, c0:c1]; disk = (yy - r) ** 2 + (xx - c) ** 2 <= rad * rad
            m = disk & (vv[r0:r1, c0:c1] > 0)
            if m.sum() < 5:
                continue
            csm = cs_now[r0:r1, c0:c1][m]; stm = stage[r0:r1, c0:c1][m]; dm = dtp[r0:r1, c0:c1][m]
            intensity = float(np.nanmean(csm)); d2p = float(np.nanmean(dm)) if np.isfinite(dm).any() else float("nan")
            series = {wl: float(np.nanmean(weekly[wl]["cs"][r0:r1, c0:c1][m])) for wl in weekly}
            timing = 1.0 if abs(d2p) <= 7 else max(0.0, 1 - (abs(d2p) - 7) / 14) if d2p > 0 else max(0.1, 1 - (abs(d2p) - 7) / 10)
            visit = 100 * (0.5 * min(intensity / 0.6, 1) + 0.25 * timing + 0.25)
            dest.append(dict(type="Feature", geometry=dict(type="Point", coordinates=[d["lon"], d["lat"]]), properties=dict(
                name=d["name"], kind=d["kind"], notes=d["notes"], radius_m=d["radius_m"], n_veg_px=int(m.sum()),
                intensity=round(intensity, 3), stage_fraction={STAGE_NAMES[i]: round(float((stm == i).mean()), 3) for i in range(1, 5)},
                days_to_peak=None if np.isnan(d2p) else round(d2p, 1), elev_m=round(float(np.nanmean(dem[r0:r1, c0:c1][m]))),
                visit_score=round(visit, 1), region=name,
                series={k_: (None if np.isnan(v) else round(v, 3)) for k_, v in series.items()})))
        json.dump(dict(type="FeatureCollection", features=dest), open(os.path.join(outdir, f"destinations_{name}.geojson"), "w"))
        print(f"[{name}] {len(dest)} destinations scored; top:", sorted([(f['properties']['name'], f['properties']['visit_score']) for f in dest], key=lambda t: -t[1])[:6])

if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "composite":
        composite(sys.argv[2], sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else "current")
    elif cmd == "phenology":
        phenology(sys.argv[2], sys.argv[3:])
    elif cmd == "region":
        region(*sys.argv[2:7])
