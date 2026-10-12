"""Build the worked example for the Method tab: image chips and chart data for a small
reference area (Washington Park Arboretum), computed from the same cached arrays the map uses.

Usage: python method_figures.py OUTDIR(docs/data/method)
"""
import os, sys, json, glob, datetime as dt
import numpy as np
from PIL import Image
from pyproj import Transformer
import config as C
from analyze import indices, colour_score, INDEX_NAMES, STAGE_NAMES, resample_to_grid

CACHE = C.SCRATCH
REF = dict(name="Washington Park Arboretum", lon=-122.2946, lat=47.6395, half_px=60)   # 120 px = 1.2 km at 10 m
STAGE_COL = {0: (20, 22, 26), 1: (26, 128, 38), 2: (242, 217, 26), 3: (250, 140, 13), 4: (120, 85, 60)}
RAMP = np.array([[26, 128, 38], [242, 217, 26], [250, 140, 13], [217, 26, 13]], float)
TEAL = np.array([[219, 230, 232], [112, 181, 194], [31, 94, 107], [8, 46, 59]], float)
BAND_NM = {"blue": 490, "green": 560, "red": 665, "rededge1": 705, "rededge2": 740, "rededge3": 783, "nir": 842, "nir08": 865, "swir16": 1610, "swir22": 2190}

def ramp(x, cols):
    x = np.clip(np.nan_to_num(x, nan=0), 0, 1) * (len(cols) - 1); i = np.clip(x.astype(int), 0, len(cols) - 2); f = (x - i)[..., None]
    return cols[i] * (1 - f) + cols[i + 1] * f

def save(name, arr, out, scale=3):
    im = Image.fromarray(np.clip(arr, 0, 255).astype(np.uint8))
    im = im.resize((im.width * scale, im.height * scale), Image.NEAREST)
    im.save(os.path.join(out, name), "WEBP", quality=88, method=4)
    return name

def rgb_chip(r, g, b):
    def st(x): return np.clip((np.nan_to_num(x, nan=0) / 0.25) ** 0.7, 0, 1) * 255
    return np.dstack([st(r), st(g), st(b)])

def main(out):
    os.makedirs(out, exist_ok=True)
    grid = json.load(open(os.path.join(CACHE, "s2_seattle", "grid.json"))); T = grid["transform"]
    tr = Transformer.from_crs(4326, C.UTM, always_xy=True); x, y = tr.transform(REF["lon"], REF["lat"])
    c = int((x - T[2]) / T[0]); r = int((y - T[5]) / T[4]); h = REF["half_px"]
    sl = (slice(r - h, r + h), slice(c - h, c + h))
    to_ll = Transformer.from_crs(C.UTM, 4326, always_xy=True)
    x0, y1 = T[2] + (c - h) * T[0], T[5] + (r - h) * T[4]; x1, y0 = x0 + 2 * h * T[0], y1 + 2 * h * T[4]
    w_, s_ = to_ll.transform(x0, y0); e_, n_ = to_ll.transform(x1, y1)
    meta = dict(ref=REF, bounds=[[s_, w_], [n_, e_]], px_m=10, size_px=2 * h, figures={})

    comp = np.load(os.path.join(CACHE, "out", "comp_seattle.npz")); ana = np.load(os.path.join(CACHE, "out", "analysis_seattle.npz"))
    anc = np.load(os.path.join(CACHE, "anc_seattle.npz")); bins = list(comp["bins"]); base = bins[0]
    veg = ana["veg"][sl]; decid = (ana["decid"][sl].astype(float) - 1) / 254; stage = ana["stage"][sl]

    # 1. weekly true colour chips (cloud-free medians) ------------------------------------------
    weeks = []
    for lab in bins:
        if f"{lab}|r" not in comp.files: continue
        rr, gg, bb = comp[f"{lab}|r"][sl], comp[f"{lab}|g"][sl], comp[f"{lab}|b"][sl]
        chip = rgb_chip(rr, gg, bb); chip[np.isnan(rr)] = (20, 22, 26)
        weeks.append(dict(label=str(lab), file=save(f"week_{lab.replace(' ', '_').replace('-', '_')}.webp", chip, out), coverage=round(float(np.isfinite(rr).mean()), 2)))
    rr, gg, bb = comp["latest|r"][sl], comp["latest|g"][sl], comp["latest|b"][sl]
    chip = rgb_chip(rr, gg, bb); chip[np.isnan(rr)] = (20, 22, 26)
    weeks.append(dict(label="Now", file=save("week_Now.webp", chip, out), coverage=round(float(np.isfinite(rr).mean()), 2)))
    meta["figures"]["weeks"] = weeks

    # 2. cloud masking on one raw pass ---------------------------------------------------------
    scenes = sorted(glob.glob(os.path.join(CACHE, "s2_seattle", "*.npz")))
    pick = None
    for f in scenes[::-1]:
        z = np.load(f); v = z["valid"][sl]
        if 0.25 < v.mean() < 0.85: pick = (f, z); break
    if pick is None: pick = (scenes[-1], np.load(scenes[-1]))
    f, z = pick; bands = z["bands"][:, sl[0], sl[1]].astype(np.float32) / C.BOA_SCALE; valid = z["valid"][sl]; scl = z["scl"][sl]
    raw = rgb_chip(bands[C.BAND_ORDER.index("red")], bands[C.BAND_ORDER.index("green")], bands[C.BAND_ORDER.index("blue")])
    masked = raw.copy(); tint = np.array([200, 40, 160], float); m = ~valid
    masked[m] = masked[m] * 0.35 + tint * 0.65
    meta["figures"]["mask"] = dict(date=os.path.basename(f)[:10], raw=save("mask_raw.webp", raw, out), masked=save("mask_applied.webp", masked, out),
                                   masked_fraction=round(float(m.mean()), 2),
                                   scl_counts={int(k): int(v) for k, v in zip(*np.unique(scl, return_counts=True))})

    # 3. spectra of a deciduous tree pixel (Sept vs Oct) and a conifer pixel ----------------------
    csn = (ana["cs_now"][sl].astype(float) - 1) / 254
    loaded = {os.path.basename(f)[:10]: np.load(f) for f in scenes}
    def spectrum(date, idx):
        z = loaded.get(date)
        if z is None or not z["valid"][sl][idx]: return None
        b = z["bands"][:, sl[0], sl[1]][:, idx[0], idx[1]].astype(float) / C.BOA_SCALE
        return {k: round(float(b[C.BAND_ORDER.index(k)]), 4) for k in C.BAND_ORDER}
    SEP = ["2026-09-16", "2026-09-08", "2026-09-18", "2026-09-06"]; OCT = ["2026-10-06", "2026-10-03", "2026-10-08", "2026-10-01"]
    def pick_pixel(cand, score):
        order = np.argsort(-np.where(cand, score, -9), axis=None)
        for flat in order[:400]:
            idx = np.unravel_index(flat, score.shape)
            if not cand[idx]: break
            d1 = next((d for d in SEP if spectrum(d, idx)), None); d2 = next((d for d in OCT if spectrum(d, idx)), None)
            if d1 and d2: return idx, d1, d2
        raise RuntimeError("no pixel with both a September and an October clear observation")
    dec_idx, sep, oct_ = pick_pixel(veg & (decid > 0.6) & anc["park_mask"][sl], csn)
    con_idx, sep_c, oct_c = pick_pixel(veg & (decid < 0.12) & anc["park_mask"][sl] & (comp[f"{base}|ndvi"][sl] > 0.6), comp[f"{base}|ndvi"][sl])
    meta["figures"]["spectra"] = dict(band_nm=BAND_NM, deciduous_sept=dict(date=sep, values=spectrum(sep, dec_idx)), deciduous_oct=dict(date=oct_, values=spectrum(oct_, dec_idx)),
                                      conifer_sept=dict(date=sep_c, values=spectrum(sep_c, con_idx)), conifer_oct=dict(date=oct_c, values=spectrum(oct_c, con_idx)),
                                      deciduous_px=[int(dec_idx[0]), int(dec_idx[1])], conifer_px=[int(con_idx[0]), int(con_idx[1])])
    # index values at those pixels, baseline vs now
    def idx_vals(lab, idx): return {k: (None if np.isnan(comp[f"{lab}|{k}"][sl][idx]) else round(float(comp[f"{lab}|{k}"][sl][idx]), 3)) for k in ("hue", "chroma", "psri", "ndvi", "cire", "ndmi")}
    last = "latest"
    meta["figures"]["index_examples"] = dict(baseline_label=base, now_label="Now (latest clear pass)", deciduous=dict(base=idx_vals(base, dec_idx), now=idx_vals(last, dec_idx)), conifer=dict(base=idx_vals(base, con_idx), now=idx_vals(last, con_idx)))

    # 4. index maps: hue, PSRI, NDVI at baseline and latest week; colour score; stage ----------------
    maps = {}
    for lab, tag in ((base, "base"), ("latest", "now")):
        hue = comp[f"{lab}|hue"][sl]; psri = comp[f"{lab}|psri"][sl]; ndvi = comp[f"{lab}|ndvi"][sl]
        maps[f"hue_{tag}"] = save(f"hue_{tag}.webp", np.where(np.isfinite(hue)[..., None], ramp(np.clip((120 - hue) / 110, 0, 1), RAMP), 24), out)   # green 120° -> red 10°
        maps[f"psri_{tag}"] = save(f"psri_{tag}.webp", np.where(np.isfinite(psri)[..., None], ramp(np.clip((psri + 0.05) / 0.35, 0, 1), TEAL), 24), out)
        maps[f"ndvi_{tag}"] = save(f"ndvi_{tag}.webp", np.where(np.isfinite(ndvi)[..., None], ramp(np.clip((ndvi - 0.1) / 0.8, 0, 1), np.array([[240, 240, 232], [120, 170, 90], [10, 70, 20]], float)), 24), out)
    cs_now = csn; maps["cs_now"] = save("cs_now.webp", np.where(veg[..., None], ramp(cs_now, RAMP), 24), out)
    st_rgb = np.zeros(stage.shape + (3,), float)
    for k, col in STAGE_COL.items(): st_rgb[stage == k] = col
    maps["stage_now"] = save("stage_now.webp", st_rgb, out)
    maps["decid"] = save("decid.webp", np.where(veg[..., None], ramp(decid, TEAL), 24), out)
    veg_rgb = np.where(veg[..., None], np.array([200, 205, 190], float), 24); maps["veg"] = save("veg.webp", veg_rgb, out)
    meta["figures"]["maps"] = maps
    vv = veg & (decid > 0.25)
    meta["figures"]["stage_fraction_ref"] = {STAGE_NAMES[i]: round(float((stage[vv] == i).mean()), 3) for i in range(1, 5)}

    # 5. phenology: baseline curve for the reference cell vs this year -------------------------------
    ph = np.load(os.path.join(CACHE, "out", "pheno.npz")); pt = ph["transform"]
    TT = np.array(T)
    curve = resample_to_grid(ph["curve"], pt, int(ph["res"]), TT, ana["veg"].shape)[:, sl[0], sl[1]]
    amp = resample_to_grid(ph["amp"], pt, int(ph["res"]), TT, ana["veg"].shape)[sl]
    peak_doy = resample_to_grid(ph["peak_doy"], pt, int(ph["res"]), TT, ana["veg"].shape)[sl]
    doys = ph["doys"].tolist()
    clim = [float(np.nanmean((curve[i] * amp)[vv])) for i in range(curve.shape[0])]
    obs = []
    for lab in bins[1:]:
        if f"{lab}|hue" not in comp.files: continue
        cs = colour_score(comp[f"{lab}|hue"][sl], comp[f"{lab}|chroma"][sl], comp[f"{lab}|psri"][sl], comp[f"{lab}|ndvi"][sl], comp[f"{base}|hue"][sl], comp[f"{base}|psri"][sl])
        rng = comp["bin_ranges"][bins.index(lab)].split("|"); mid = dt.date(2026, int(rng[0][:2]), int(rng[0][3:])) + dt.timedelta(days=3)
        obs.append(dict(label=str(lab), doy=mid.timetuple().tm_yday, value=round(float(np.nanmean(cs[vv])), 3), coverage=round(float(np.isfinite(cs[vv]).mean()), 2)))
    fc = []
    for rec in ana["fc_bins"]:
        lab, a0, a1 = str(rec).split("|"); mid = dt.date.fromisoformat(a0) + dt.timedelta(days=3)
        v = (ana[f"fc|{lab}|cs"][sl].astype(float) - 1) / 254; v[ana[f"fc|{lab}|cs"][sl] == 0] = np.nan
        fc.append(dict(label=lab, doy=mid.timetuple().tm_yday, value=round(float(np.nanmean(v[vv])), 3)))
    meta["figures"]["phenology"] = dict(doys=doys, climatology=[round(v, 3) for v in clim], observed=obs, forecast=fc,
                                        peak_doy_clim=round(float(np.nanmean(peak_doy[vv])), 1), amp_clim=round(float(np.nanmean(amp[vv])), 3),
                                        year_doy_origin="2026-01-01")
    # 6. weather: SeaTac nightly lows, this year vs 2025 same dates, plus GFS forecast ----------------
    w = json.load(open(os.path.join(CACHE, "weather", "weather.json")))
    st26 = w["stations"]["2026"]["USW00024233"]; st25 = w["stations"]["2025"]["USW00024233"]
    lows26 = [(d[4:], v[0]) for d, v in sorted(st26["series"].items()) if d >= "20260915" and v[0] is not None]
    lows25 = [(d[4:], v[0]) for d, v in sorted(st25["series"].items()) if "0915" <= d[4:] <= "1115" and v[0] is not None]
    gfs = w["gfs"]["Seattle"]
    meta["figures"]["weather"] = dict(station=st26["name"], lows_2026=lows26, lows_2025=lows25, forecast=[(d["date"][5:].replace("-", ""), d["tmin"]) for d in gfs], gfs_cycle=w["gfs_cycle"],
                                      nights_le_7_2026=st26["nights_le_7"], nights_le_7_2025_same_dates=sum(1 for _, v in lows25 if v <= 7 and _ <= lows26[-1][0]))
    # 7. hotspot aggregation: 100 m cells of colour score over the chip, with the threshold used --------
    k = 10; Hc, Wc = (2 * h) // k, (2 * h) // k
    cs_cells = np.zeros((Hc, Wc)); frac = np.zeros((Hc, Wc))
    for i in range(Hc):
        for j in range(Wc):
            blk = vv[i * k:(i + 1) * k, j * k:(j + 1) * k]; v = cs_now[i * k:(i + 1) * k, j * k:(j + 1) * k]
            frac[i, j] = blk.mean(); cs_cells[i, j] = float(np.nanmean(v[blk])) if blk.any() else 0
    score = cs_cells * np.clip(frac / 0.25, 0, 1)
    cell_rgb = ramp(score, RAMP); cell_rgb[frac < 0.2] = 40
    meta["figures"]["cells"] = dict(file=save("cells.webp", np.kron(cell_rgb, np.ones((k, k, 1))), out, scale=3), cell_m=100,
                                    scores=[[round(float(s), 2) for s in row] for row in score], fractions=[[round(float(s), 2) for s in row] for row in frac],
                                    threshold_note="top tenth of cells county-wide, with at least 20% tree cover")
    json.dump(meta, open(os.path.join(out, "method.json"), "w"), indent=1)
    print("ref pixel deciduous", dec_idx, "conifer", con_idx, "| masked scene", meta["figures"]["mask"]["date"], meta["figures"]["mask"]["masked_fraction"])
    print("stage fractions in reference area", meta["figures"]["stage_fraction_ref"]); print("phenology clim peak DOY", meta["figures"]["phenology"]["peak_doy_clim"], "obs", obs, "fc", fc)
    print("spectra dec sept/oct", sep, oct_, meta["figures"]["spectra"]["deciduous_sept"]["values"]["red"], meta["figures"]["spectra"]["deciduous_oct"]["values"]["red"])

if __name__ == "__main__":
    main(sys.argv[1])
