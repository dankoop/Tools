# Ground-truth check: southern Shoreline, drone ortho of 2026-10-11 vs the map's Now layer (Sentinel-2, 2026-10-08)

Registration: best match at zero shift (edge correlation 0.48, brightness 0.73), so alignment is within one 10 m cell. Oct 8 pass is not hazy (blue 0.054 vs 0.041-0.078 in September passes).

## What the predicted cells actually contain (labelled by eye from the 2 cm ortho)

| predicted | peak_decid | turning_decid | green | conifer | lawn | nonveg | n |
|---|---|---|---|---|---|---|---|
| peak | 3 | 2 | 4 | 8 | 0 | 7 | 24 |
| turning | 0 | 2 | 6 | 6 | 1 | 9 | 24 |
| green | 0 | 1 | 4 | 4 | 1 | 14 | 24 |
| evergreen_hist | 0 | 0 | 2 | 2 | 3 | 9 | 16 |
| nonveg | 0 | 0 | 0 | 0 | 0 | 8 | 8 |

Precision of 'peak' cells (any turning or peak deciduous tree inside): 5/24. 'Turning': 2/24. 31 of 48 cells predicted turning or peak hold no deciduous tree at all.

Recall sample (32 cells chosen where the drone image shows warm colour, verified by eye): 7 held real colour; the map called 1 peak, 4 green with deciduous likelihood near 0, 2 were masked as non-vegetation.

Mechanism: conifer canopy is dark (reflectance sum 0.07-0.09 vs 0.15 for turned trees). Hue computed from such small values is unstable, and the lower October sun shifted conifer hue from ~120 to ~95 degrees with no real colour change, which the score reads as turning. Built, mixed cells pass the September NDVI > 0.4 vegetation test and then produce spurious hue shifts.

Tested fix on the labelled cells: gating the hue term by baseline brightness cuts false 'peak' cells from 18 to 1 of 87 at the 0.45 threshold while keeping 2 of 4 true peaks; AUC 0.71 -> 0.77. Not yet applied to the map.
