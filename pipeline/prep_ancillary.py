"""Prepare ancillary layers on the analysis grids and summarise weather.

Outputs (in CACHE):
  anc_<region>.npz      landcover (ESA WorldCover class), dem (m), kc_mask, seattle_mask, nhood_id
  weather/weather.json  station series + metrics, HRRR recent-day summaries, GFS 7-day forecast
  weather/hrrr_daily.npz  per-day Tmin/Tmax/gust grids on a 1 km county grid (for the analysis step)

Usage: python prep_ancillary.py
"""
import os, json, glob, csv, datetime as dt, collections
import numpy as np, rasterio
from rasterio.warp import reproject, Resampling
from rasterio import features
from rasterio.transform import Affine
from shapely.geometry import shape, mapping
from shapely.ops import unary_union
from pyproj import Transformer
import config as C

CACHE = C.SCRATCH
W = os.path.join(CACHE, "weather")
B = os.path.join(CACHE, "bounds")
A = os.path.join(CACHE, "anc")

def load_grid(region_dir):
    g = json.load(open(os.path.join(region_dir, "grid.json")))
    return Affine(*g["transform"]), g["width"], g["height"], g["crs"]

def warp_to_grid(src_path_or_list, transform, width, height, crs, resampling, dtype, nodata=0):
    out = np.full((height, width), nodata, dtype)
    paths = src_path_or_list if isinstance(src_path_or_list, list) else [src_path_or_list]
    for p in paths:
        with rasterio.open(p) as src:
            tmp = np.full((height, width), nodata, dtype)
            reproject(rasterio.band(src, 1), tmp, dst_transform=transform, dst_crs=crs,
                      resampling=resampling, dst_nodata=nodata, src_nodata=src.nodata)
            m = tmp != nodata
            out[m] = tmp[m]
    return out

def geoms_utm(path, to_crs):
    tr = Transformer.from_crs(4326, to_crs, always_xy=True)
    import shapely
    fc = json.load(open(path))
    out = []
    for f in fc["features"]:
        g = shape(f["geometry"])
        out.append((shapely.transform(g, lambda xy: np.column_stack(tr.transform(xy[:, 0], xy[:, 1]))), f.get("properties", {})))
    return out

def prep_region(region_dir, name):
    transform, width, height, crs = load_grid(region_dir)
    print(f"[{name}] grid {width}x{height}")
    lc = warp_to_grid(os.path.join(A, "worldcover_N45W123.tif"), transform, width, height, crs, Resampling.nearest, np.uint8)
    dem = warp_to_grid(sorted(glob.glob(os.path.join(A, "dem_*.tif"))), transform, width, height, crs, Resampling.bilinear, np.float32, nodata=-9999)
    kc = unary_union([g for g, _ in geoms_utm(os.path.join(B, "king_county.geojson"), crs)])
    sea = unary_union([g for g, _ in geoms_utm(os.path.join(B, "city-limits.geojson"), crs)])
    kc_mask = features.rasterize([(mapping(kc), 1)], out_shape=(height, width), transform=transform, dtype=np.uint8).astype(bool)
    sea_mask = features.rasterize([(mapping(sea), 1)], out_shape=(height, width), transform=transform, dtype=np.uint8).astype(bool)
    nh = geoms_utm(os.path.join(B, "neighborhoods.geojson"), crs)
    nh_names = sorted({p.get("name") or p.get("nhood") for _, p in nh})
    nh_id = features.rasterize([(mapping(g), nh_names.index(p.get("name") or p.get("nhood")) + 1) for g, p in nh],
                               out_shape=(height, width), transform=transform, dtype=np.uint16)
    parks = geoms_utm(os.path.join(B, "seattle-parks-osm.geojson"), crs)
    park_mask = features.rasterize([(mapping(g), 1) for g, _ in parks], out_shape=(height, width), transform=transform, dtype=np.uint8).astype(bool)
    np.savez_compressed(os.path.join(CACHE, f"anc_{name}.npz"), landcover=lc, dem=dem, kc_mask=kc_mask,
                        seattle_mask=sea_mask, nhood_id=nh_id, park_mask=park_mask, nhood_names=np.array(nh_names))
    print(f"[{name}] landcover classes {dict(zip(*np.unique(lc[kc_mask], return_counts=True)))}")
    print(f"[{name}] dem range {np.nanmin(dem[kc_mask]):.0f}-{np.nanmax(dem[kc_mask]):.0f} m; parks px {park_mask.sum()}")

# ---------------------------------------------------------------- weather
def parse_stations():
    meta = {}
    for l in open(os.path.join(W, "ghcnd-stations.txt")):
        meta[l[:11]] = dict(lat=float(l[12:20]), lon=float(l[21:30]), elev=float(l[31:37]), name=l[41:71].strip().title())
    return meta

def read_ghcn(years):
    data = collections.defaultdict(lambda: collections.defaultdict(dict))  # sid -> element -> date -> value
    for y in years:
        p = os.path.join(W, f"ghcn_kc_{y}.csv")
        if not os.path.exists(p):
            continue
        for r in csv.reader(open(p)):
            if len(r) >= 6 and r[2] in ("TMIN", "TMAX", "PRCP", "AWND", "WSF2", "WSF5") and r[5] == "":
                data[r[0]][r[2]][r[1]] = int(r[3])
    return data

def station_summary(data, meta, year, since="0915", until="1231"):
    """Per-station metrics for an autumn window of a given year."""
    out = {}
    for sid, el in data.items():
        tmin = {d: v / 10 for d, v in el.get("TMIN", {}).items() if d[:4] == str(year) and since <= d[4:] <= until}
        if len(tmin) < 10:
            continue
        tmax = {d: v / 10 for d, v in el.get("TMAX", {}).items() if d[:4] == str(year)}
        prcp = {d: v / 10 for d, v in el.get("PRCP", {}).items() if d[:4] == str(year)}
        wsf = {d: v / 10 for d, v in el.get("WSF2", {}).items() if d[:4] == str(year)}
        out[sid] = dict(
            **meta.get(sid, {}), last_obs=max(tmin), n=len(tmin),
            nights_le_7=sum(v <= 7 for v in tmin.values()), nights_le_2=sum(v <= 2 for v in tmin.values()),
            nights_le_0=sum(v <= 0 for v in tmin.values()), min_tmin=min(tmin.values()),
            mean_tmin=round(float(np.mean(list(tmin.values()))), 1),
            first_frost=min([d for d, v in tmin.items() if v <= 0], default=None),
            series={d: [tmin.get(d), tmax.get(d), prcp.get(d), wsf.get(d)] for d in sorted(set(tmin) | set(tmax))
                    if d[4:] >= "0901"},
        )
    return out

def grib_daily(outdir_1km):
    """Aggregate HRRR analyses to local-day Tmin/Tmax/gust on a 1 km county grid."""
    files = sorted(glob.glob(os.path.join(W, "hrrr", "tmp2m_*.grib2")))
    gfiles = sorted(glob.glob(os.path.join(W, "hrrr", "gust_*.grib2")))
    res = 1000
    tr = Transformer.from_crs(4326, C.UTM, always_xy=True)
    lon0, lat0, lon1, lat1 = C.REGIONS["county"]["bbox"]
    x0, y0 = tr.transform(lon0, lat0); x1, y1 = tr.transform(lon1, lat1)
    x0, y0 = np.floor(x0 / res) * res, np.floor(y0 / res) * res; x1, y1 = np.ceil(x1 / res) * res, np.ceil(y1 / res) * res
    transform = Affine(res, 0, x0, 0, -res, y1); width, height = int((x1 - x0) / res), int((y1 - y0) / res)

    def warp(p):
        return warp_to_grid(p, transform, width, height, C.UTM, Resampling.bilinear, np.float32, nodata=-9999)

    # local date assignment: the night of local date D is UTC D 04Z..16Z; afternoon of D is D 22Z and D+1 00Z
    tmin, tmax, gust = {}, {}, {}
    for p in files:
        ds, hz = os.path.basename(p)[6:14], int(os.path.basename(p)[15:17])
        a = warp(p)
        if 4 <= hz <= 16:
            tmin.setdefault(ds, []).append(a)
        elif hz == 22:
            tmax.setdefault(ds, []).append(a)
        elif hz == 0:
            prev = (dt.datetime.strptime(ds, "%Y%m%d") - dt.timedelta(days=1)).strftime("%Y%m%d")
            tmax.setdefault(prev, []).append(a)
    for p in gfiles:
        ds, hz = os.path.basename(p)[5:13], int(os.path.basename(p)[14:16])
        loc = (dt.datetime.strptime(ds, "%Y%m%d") + dt.timedelta(hours=hz - 7)).strftime("%Y%m%d")
        gust.setdefault(loc, []).append(warp(p))
    days = sorted(d for d in tmin if len(tmin[d]) >= 4)
    TMIN = np.stack([np.min(tmin[d], axis=0) for d in days])
    TMAX = np.stack([np.max(tmax[d], axis=0) if d in tmax else np.full((height, width), np.nan, np.float32) for d in days])
    GUST = np.stack([np.max(gust[d], axis=0) if d in gust else np.full((height, width), np.nan, np.float32) for d in days])
    np.savez_compressed(os.path.join(W, "hrrr_daily.npz"), days=np.array(days), tmin=TMIN, tmax=TMAX, gust=GUST,
                        transform=np.array(list(transform)[:6]), width=width, height=height)
    return days, TMIN, TMAX, GUST, transform

def gfs_forecast(points):
    """Daily (local) forecast at points from the 3-hourly GFS messages."""
    meta = json.load(open(os.path.join(W, "weather_meta.json")))
    cyc = dt.datetime.strptime(meta["gfs_cycle"], "%Y%m%dT%HZ")
    rows = collections.defaultdict(lambda: collections.defaultdict(list))
    prev_apcp = {}
    for f in range(3, 169, 3):
        t_local = cyc + dt.timedelta(hours=f - 7)
        day = t_local.strftime("%Y-%m-%d")
        vals = {}
        for var in ("tmp2m", "gust", "apcp"):
            p = os.path.join(W, "gfs", f"{var}_f{f:03d}.grib2")
            if not os.path.exists(p):
                continue
            with rasterio.open(p) as ds:
                a = ds.read(1)
                for name, (lon, lat) in points.items():
                    r, c = ds.index(lon, lat)
                    vals[(var, name)] = float(a[r, c])
        for name in points:
            if ("tmp2m", name) in vals:
                rows[name][day].append(("t", vals[("tmp2m", name)]))
            if ("gust", name) in vals:
                rows[name][day].append(("g", vals[("gust", name)]))
            if ("apcp", name) in vals:
                v = vals[("apcp", name)]
                inc = v if f % 6 == 3 else v - prev_apcp.get(name, 0.0)   # GFS APCP resets every 6 h
                prev_apcp[name] = v if f % 6 == 3 else 0.0
                rows[name][day].append(("p", max(inc, 0.0)))
    out = {}
    for name, days in rows.items():
        out[name] = []
        for day in sorted(days):
            t = [v for k, v in days[day] if k == "t"]; g = [v for k, v in days[day] if k == "g"]; p = [v for k, v in days[day] if k == "p"]
            if len(t) < 4:
                continue  # partial day
            out[name].append(dict(date=day, tmin=round(min(t), 1), tmax=round(max(t), 1),
                                  gust_max_ms=round(max(g), 1) if g else None, precip_mm=round(sum(p), 1)))
    return out, meta["gfs_cycle"]

POINTS = {"Seattle": (-122.33, 47.62), "Bellevue / Eastside": (-122.15, 47.6), "Snoqualmie Valley (North Bend)": (-121.78, 47.5),
          "Snoqualmie Pass": (-121.42, 47.42), "Kent / Green River Valley": (-122.23, 47.38)}

def main():
    import sys
    for name, d in ([] if "--weather-only" in sys.argv else [("county", "s2_county"), ("seattle", "s2_seattle"), ("baseline", "s2_baseline_2025"), ("shoreline", "s2_shoreline")]):
        if os.path.exists(os.path.join(CACHE, d, "grid.json")):
            prep_region(os.path.join(CACHE, d), name)
    meta = parse_stations()
    data = read_ghcn([2023, 2024, 2025, 2026])
    stations = {y: station_summary(data, meta, y, until="1115" if y < 2026 else "1231") for y in (2023, 2024, 2025, 2026)}
    days, TMIN, TMAX, GUST, transform = grib_daily(W)
    with np.errstate(all="ignore"):
        hrrr = dict(days=[f"{d[:4]}-{d[4:6]}-{d[6:]}" for d in days],
                    county_mean_tmin=[round(float(np.nanmean(a)), 1) for a in TMIN],
                    county_mean_tmax=[None if np.all(np.isnan(a)) else round(float(np.nanmean(a)), 1) for a in TMAX],
                    county_max_gust_ms=[round(float(np.nanmax(a)), 1) for a in GUST],
                    lowland_mean_tmin=None)
    gfs, cycle = gfs_forecast(POINTS)
    json.dump(dict(stations=stations, hrrr=hrrr, gfs=gfs, gfs_cycle=cycle, points=POINTS),
              open(os.path.join(W, "weather.json"), "w"), indent=1)
    print("HRRR days", days[0], "->", days[-1], "county mean Tmin last 5:", hrrr["county_mean_tmin"][-5:])
    print("GFS", cycle, json.dumps(gfs["Seattle"][:3]))
    print("stations 2026:", {k: (v["name"], v["last_obs"], v["nights_le_2"], v["min_tmin"]) for k, v in list(stations[2026].items())[:6]})

if __name__ == "__main__":
    main()
