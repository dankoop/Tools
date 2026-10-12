# Seattle Fall Color Sentinel

A fall-colour map for Seattle and King County built from Sentinel-2 imagery, prior-year tree
phenology and NOAA weather. The web map (`docs/index.html`) shows where colour is turning now,
how it has spread week by week since early September, where peak is predicted, and which
hotspots and named destinations are worth visiting this week. A curated, toggleable layer of
dated local reports sits on top.

## What's in the repo

```
pipeline/
  config.py          regions, grids, band list
  inventory.py       list Sentinel-2 L2A scenes (AWS sentinel-cogs bucket) with STAC metadata
  ingest_s2.py       windowed COG reads -> cloud-masked band stacks on a fixed UTM grid
  fetch_weather.py   HRRR analyses + GFS 7-day forecast by GRIB byte-range
  prep_ancillary.py  WorldCover land cover, Copernicus DEM, boundaries, GHCN stations, weather summary
  analyze.py         indices, weekly composites, colour score & stages, phenology, clusters, scores
  render.py          Web-Mercator PNG/JPEG overlays + legends
  export_site.py     assembles docs/ (manifest, GeoJSON, overlays, page)
  data/places.json   named destinations and reference places
  data/reports.json  curated, dated reports with sources
  templates/         the single-page web map
docs/                the static site (GitHub Pages: Settings -> Pages -> main /docs)
```

## Method in brief

* **Imagery.** Every Sentinel-2 L2A scene over MGRS tiles 10TET and 10TFT from Sept 1 onward,
  read directly from `s3://sentinel-cogs` (no API key). Seattle at 10 m, King County at 20 m.
  Clouds, cirrus, shadow and snow masked per pixel from the SCL band plus blue/SWIR brightness
  tests, dilated 2 px, then a temporal test: a pixel whose blue reflectance sits more than 0.04
  above its own early-September median is residual cloud or haze and is dropped too. Weekly layers are per-pixel medians of clear observations; "Now" is the
  latest clear observation, back-filled from the last two weeks.
* **Spectral colour.** Hue and chroma from red/green/blue reflectance; PSRI
  `(red - blue) / red-edge-2` for senescence; NDVI and red-edge chlorophyll index for greenness
  loss; NDMI for leaf-off. The colour score blends hue shift and PSRI rise against each pixel's
  own early-September baseline, scaled down once NDVI shows the leaves are gone. Stages
  (green / turning / peak / past peak) are thresholds on that score and on the NDVI drop.
* **Deciduous test.** A pixel can be called turning or peak only if it is deciduous: its winter
  (Jan-Feb) NDVI must fall below about 0.98 of its summer (Jun-Aug) NDVI, measured from the same
  archive on the same grid. Conifers and winter-green lawns stay near or above 1; pixels that are
  not green canopy in midsummer score zero. Where no clear winter or summer pass exists the
  2023-25 autumn history is the fallback. Against the Shoreline drone labels this lifted
  turning/peak precision from 32% to 50% and cut conifer false positives from 32 to 4 of the
  flagged cells, at the cost of half the recall (partial crowns in winter-green pixels).
* **Phenology.** Autumn 2023-2025 processed the same way at 80 m gives each pixel a typical
  colour curve, peak week and amplitude (small amplitude = evergreen). This year's progress is
  matched to that curve to get a phase offset and a predicted peak date, pulled earlier when
  cold nights are forecast (more at elevation).
* **Weather.** GHCN daily stations (airports, SNOTEL) for nightly lows and cold-night counts
  versus prior years; HRRR 3 km analyses fill the days after the station feed; the latest GFS
  run gives the 7-day outlook and the leaf-drop risk (gusts, heavy rain).
* **Hotspots and scores.** Tree pixels aggregated to 100 m (Seattle) / 200 m (county) cells;
  the top decile by colour score grouped into contiguous hotspots. Visit score = 45% colour
  intensity, 20% timing against predicted peak, 20% access (parks, nearby neighbourhoods),
  15% size. Named destinations are scored on the pixels inside their footprint.

## Re-running

Needs Python 3.11+ and network access to AWS S3 (public buckets) and
raw.githubusercontent.com (boundaries). `pip install -r pipeline/requirements.txt`.

```bash
export FALLCOLOR_CACHE=/path/to/big/cache     # ~10-15 GB for a full season
cd pipeline
python inventory.py $FALLCOLOR_CACHE/scenes_2026.json 2026 9 10 11
python ingest_s2.py county  $FALLCOLOR_CACHE/scenes_2026.json $FALLCOLOR_CACHE/s2_county  --min-date 2026-09-01
python ingest_s2.py seattle $FALLCOLOR_CACHE/scenes_2026.json $FALLCOLOR_CACHE/s2_seattle --min-date 2026-09-01
# baseline years (once)
for y in 2023 2024 2025; do
  python inventory.py $FALLCOLOR_CACHE/scenes_$y.json $y 9 10 11
  python ingest_s2.py county $FALLCOLOR_CACHE/scenes_$y.json $FALLCOLOR_CACHE/s2_baseline_$y --res 80 --min-date $y-09-01 --max-date $y-11-30
  python analyze.py composite $FALLCOLOR_CACHE/s2_baseline_$y $FALLCOLOR_CACHE/out/comp_base_$y.npz baseline
done
python analyze.py phenology $FALLCOLOR_CACHE/out/pheno.npz $FALLCOLOR_CACHE/out/comp_base_202{3,4,5}.npz
# weather + ancillary
python fetch_weather.py $FALLCOLOR_CACHE/weather 20260915 $(date +%Y%m%d)
#   (GHCN: stream csv/by_year/<year>.csv from the noaa-ghcn-pds bucket filtered to King County station ids
#    into $FALLCOLOR_CACHE/weather/ghcn_kc_<year>.csv; see prep_ancillary.py)
python prep_ancillary.py
# winter and summer passes for the deciduous test (once a year)
python inventory.py $FALLCOLOR_CACHE/scenes_winter.json 2026 1 2      # plus the previous January
python inventory.py $FALLCOLOR_CACHE/scenes_summer.json 2026 6 7 8
for r in seattle county; do
  python ingest_s2.py $r $FALLCOLOR_CACHE/scenes_winter.json $FALLCOLOR_CACHE/s2_winter_$r --min-date 2025-01-01 --max-cloud 20
  python ingest_s2.py $r $FALLCOLOR_CACHE/scenes_summer.json $FALLCOLOR_CACHE/s2_summer_$r --min-date 2026-06-01 --max-cloud 10
  python analyze.py season $r $FALLCOLOR_CACHE/out
done
# analysis, rendering, site
for r in seattle county; do
  python analyze.py composite $FALLCOLOR_CACHE/s2_$r $FALLCOLOR_CACHE/out/comp_$r.npz current
  python analyze.py region $r $FALLCOLOR_CACHE/s2_$r $FALLCOLOR_CACHE/out/comp_$r.npz $FALLCOLOR_CACHE/out/pheno.npz $FALLCOLOR_CACHE/out
  python render.py $r $FALLCOLOR_CACHE/out/analysis_$r.npz $FALLCOLOR_CACHE/out/render
done
# sharp 10 m colour bases: latest clear pixel (the greyscale default) and the cloud-free Sep 16 "before" image
python base_imagery.py $FALLCOLOR_CACHE/scenes_2026.json $FALLCOLOR_CACHE/out --min-date 2026-09-15 --max-cloud 60
python base_imagery.py $FALLCOLOR_CACHE/scenes_2026.json $FALLCOLOR_CACHE/out --min-date 2026-09-01 --dates 2026-09-16,2026-09-08 --name base_sep16 --label "Sep 16 (before)"
python export_site.py $FALLCOLOR_CACHE/out ../docs
```

Ingest is incremental (dates already on disk are skipped), so a refresh re-reads only new scenes.

## Data sources

Copernicus Sentinel-2 L2A (ESA), COGs by Element 84 on AWS Open Data · ESA WorldCover 2021 ·
Copernicus DEM GLO-30 · NOAA GHCN-Daily, HRRR and GFS on AWS Open Data · Seattle boundaries,
neighbourhoods and parks from seattleio/seattle-boundaries-data (OSM / City of Seattle) ·
King County outline from US Census via plotly datasets. Reports are quoted from and linked to
their original publishers.
