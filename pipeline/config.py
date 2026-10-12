"""Shared configuration for the Seattle / King County fall-color pipeline."""
import os

SCRATCH = os.environ.get("FALLCOLOR_CACHE", os.path.join(os.path.dirname(__file__), "..", "cache"))
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DOCS = os.path.join(REPO, "docs")

BUCKET = "https://sentinel-cogs.s3.us-west-2.amazonaws.com"
UTM = "EPSG:32610"

# Bands read for every scene (asset key -> native resolution, m)
BANDS = {"blue": 10, "green": 10, "red": 10, "rededge1": 20, "rededge2": 20, "rededge3": 20,
         "nir": 10, "nir08": 20, "swir16": 20, "swir22": 20}
BAND_ORDER = list(BANDS)
# Reflectance = (DN - BOA_OFFSET) / BOA_SCALE
# Element84 sentinel-cogs items carry earthsearch:boa_offset_applied=True, i.e. the stored DN is
# already reflectance*10000 (verified against the data: red DN ~100-1000 over forest).
BOA_OFFSET, BOA_SCALE = 0, 10000.0

# Analysis regions: name -> (lon/lat bbox, resolution m).  Bboxes are snapped to the grid.
REGIONS = {
    # King County bounding box with ~1 km margin
    "county": dict(bbox=(-122.55, 47.07, -121.05, 47.79), res=20),
    # Seattle city limits with ~500 m margin
    "seattle": dict(bbox=(-122.445, 47.48, -122.22, 47.745), res=10),
    # validation window: southern Shoreline drone orthoimage footprint with a 300 m margin
    "shoreline": dict(bbox=(-122.358, 47.7405, -122.336, 47.7565), res=10),
}

# SCL classes kept as potentially-valid land observations
SCL_KEEP = {4, 5, 7}          # vegetation, not-vegetated, unclassified
SCL_WATER = 6
