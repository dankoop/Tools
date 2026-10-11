"""Fetch NOAA HRRR analyses (recent 2 m temperature / surface gust) and the latest GFS
7-day forecast for the King County area via byte-range reads of the public S3 GRIB2 files.

Each GRIB2 message is saved as its own small .grib2 file (readable by rasterio/GDAL).
Usage: python fetch_weather.py OUTDIR START_DATE END_DATE   (dates YYYYMMDD, inclusive)
"""
import os, sys, re, datetime as dt, concurrent.futures as cf, urllib.request, json

HRRR = "https://noaa-hrrr-bdp-pds.s3.amazonaws.com"
GFS = "https://noaa-gfs-bdp-pds.s3.amazonaws.com"
# HRRR analysis hours (UTC) to sample: overnight/morning for Tmin, afternoon for Tmax
HRRR_TMP_HOURS = [4, 6, 8, 10, 12, 14, 16, 22, 0]
HRRR_GUST_HOURS = [0, 3, 6, 9, 12, 15, 18, 21]

def get(url, timeout=120):
    return urllib.request.urlopen(url, timeout=timeout).read()

def idx_ranges(idx_text):
    """Parse a .idx file -> list of (start, end_or_None, description)."""
    lines = [l for l in idx_text.splitlines() if l.strip()]
    out = []
    for i, l in enumerate(lines):
        parts = l.split(":")
        start = int(parts[1])
        end = int(lines[i + 1].split(":")[1]) - 1 if i + 1 < len(lines) else None
        out.append((start, end, ":" + ":".join(parts[3:])))
    return out

def fetch_message(grib_url, want, out_path):
    """Download the GRIB message whose idx description contains all strings in `want`."""
    if os.path.exists(out_path) and os.path.getsize(out_path) > 1000:
        return out_path
    try:
        idx = get(grib_url + ".idx").decode()
    except Exception as e:
        return f"MISSING {grib_url}: {e}"
    for start, end, desc in idx_ranges(idx):
        if all(w in desc for w in want):
            req = urllib.request.Request(grib_url, headers={"Range": f"bytes={start}-{end if end else ''}"})
            data = urllib.request.urlopen(req, timeout=300).read()
            with open(out_path, "wb") as f:
                f.write(data)
            return out_path
    return f"NOFIELD {grib_url} {want}"

def hrrr_jobs(outdir, d0, d1):
    jobs = []
    d = d0
    while d <= d1:
        ds = d.strftime("%Y%m%d")
        for h in HRRR_TMP_HOURS:
            url = f"{HRRR}/hrrr.{ds}/conus/hrrr.t{h:02d}z.wrfsfcf00.grib2"
            jobs.append((url, [":TMP:2 m above ground:"], f"{outdir}/hrrr/tmp2m_{ds}_{h:02d}z.grib2"))
        for h in HRRR_GUST_HOURS:
            url = f"{HRRR}/hrrr.{ds}/conus/hrrr.t{h:02d}z.wrfsfcf00.grib2"
            jobs.append((url, [":GUST:surface:"], f"{outdir}/hrrr/gust_{ds}_{h:02d}z.grib2"))
        d += dt.timedelta(days=1)
    return jobs

def latest_gfs_cycle():
    now = dt.datetime.utcnow()
    for back in range(0, 8):
        t = now - dt.timedelta(hours=6 * back)
        t = t.replace(hour=(t.hour // 6) * 6)
        ds, hh = t.strftime("%Y%m%d"), t.strftime("%H")
        url = f"{GFS}/gfs.{ds}/{hh}/atmos/gfs.t{hh}z.pgrb2.0p25.f168.idx"
        try:
            get(url, timeout=30)
            return ds, hh
        except Exception:
            continue
    raise RuntimeError("no complete GFS cycle found")

def gfs_jobs(outdir, ds, hh):
    jobs = []
    for f in range(3, 169, 3):
        url = f"{GFS}/gfs.{ds}/{hh}/atmos/gfs.t{hh}z.pgrb2.0p25.f{f:03d}"
        jobs.append((url, [":TMP:2 m above ground:"], f"{outdir}/gfs/tmp2m_f{f:03d}.grib2"))
        jobs.append((url, [":GUST:surface:"], f"{outdir}/gfs/gust_f{f:03d}.grib2"))
        jobs.append((url, [":APCP:surface:"], f"{outdir}/gfs/apcp_f{f:03d}.grib2"))
    return jobs

def main(outdir, start, end):
    os.makedirs(f"{outdir}/hrrr", exist_ok=True); os.makedirs(f"{outdir}/gfs", exist_ok=True)
    d0 = dt.datetime.strptime(start, "%Y%m%d"); d1 = dt.datetime.strptime(end, "%Y%m%d")
    ds, hh = latest_gfs_cycle()
    json.dump({"gfs_cycle": f"{ds}T{hh}Z", "hrrr_start": start, "hrrr_end": end},
              open(f"{outdir}/weather_meta.json", "w"))
    jobs = hrrr_jobs(outdir, d0, d1) + gfs_jobs(outdir, ds, hh)
    print(f"{len(jobs)} GRIB messages to fetch; GFS cycle {ds} {hh}Z")
    bad = 0
    with cf.ThreadPoolExecutor(12) as ex:
        for r in ex.map(lambda j: fetch_message(*j), jobs):
            if not r.endswith(".grib2"):
                bad += 1; print(r)
    print("done; problems:", bad)

if __name__ == "__main__":
    main(*sys.argv[1:4])
