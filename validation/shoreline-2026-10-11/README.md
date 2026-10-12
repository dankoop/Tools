# Ground-truth check: southern Shoreline, drone ortho of 2026-10-11 vs the map's Now layer (Sentinel-2, 2026-10-08)

Registration: best match at zero shift (edge correlation 0.48, brightness 0.73), so alignment is within one 10 m cell. Oct 8 pass is not hazy (blue 0.054 vs 0.041-0.078 in September passes).

## Full check: 289 cells labelled by the image owner (Oct 12), map's Now layer from Oct 8

Cells were sampled across the map's predicted classes plus 40 cells that looked warm-coloured in the drone image, shown blind (prediction hidden), labelled by what fills the 10 m square and how much of it a deciduous crown covers.

| map prediction | peak | turning | green deciduous | conifer | lawn | built | n |
|---|---|---|---|---|---|---|---|
| peak | 6 | 12 | 5 | 22 | 5 | 10 | 60 |
| turning | 2 | 10 | 8 | 27 | 8 | 5 | 60 |
| green | 0 | 13 | 2 | 27 | 8 | 14 | 64 |
| masked | 3 | 18 | 4 | 21 | 28 | 31 | 105 |

**Precision.** Of 60 cells the map calls peak, 18 hold a turning or peak deciduous tree (30%); of 60 called turning, 12 do (20%). 101 of the 120 cells called turning or peak hold no colour; most are conifer canopy (49) or lawn, roofs and pavement (28).

**Recall.** Of the 64 cells you labelled as holding a turning or peak tree, the map calls 18 peak, 12 turning, 13 green and excludes 21 through its vegetation and deciduous masks. Only 18 of the 64 have a crown covering more than half the cell; 43 are partial crowns.

**Can tuning fix it?** No, not with these inputs. On a geographic split (alternating 200 m blocks), every single satellite feature separates colour cells from the rest with AUC 0.52 to 0.58, a 16-feature logistic model reaches 0.55 held-out, and the best hand rule 0.65. Restricting to crowns that fill most of the cell helps only a little (AUC 0.72). For the 11 cells labelled peak, the best rule (warm hue, hue shift, and a brightness gate) reaches AUC 0.74 and 14% precision at 73% recall. The brightness gate does remove most conifer false positives but loses real cells too; the net change is within noise, so the map's rules were left as they are.

**Reading.** In a mixed suburb a 10 m pixel is a blend of a crown, its neighbours, a roof and shadow, and the three-day gap to the drone pass does not explain the result. The layer should be read as indicative at pixel level in residential areas, and the ranked places and regional summaries, which aggregate many pixels, are the trustworthy products. Three things would change this: a tree-canopy and species layer (King County LiDAR canopy, city tree inventories) to replace the spectral masks; aggregation to 30 to 50 m before classification; or imagery that resolves crowns (NAIP at 60 cm, PlanetScope at 3 m). A second drone flight over a park or forested area would show whether the layer performs where crowns fill cells.

Files: `labelled_cells_290.json` (labels joined to every satellite value), chip sheets and the overlay above, and the labelling page under `/labeler` in the repository.

## Re-scored on Oct 12 against the rebuilt layer (Sentinel-2 pass of 2026-10-11, same day as the drone flight)

| layer | flagged turning or peak | hold a turning or peak tree | conifer among flagged | lawn or built among flagged | colour cells found (of 64) | full-crown colour cells found (of 18) |
|---|---|---|---|---|---|---|
| Oct 11 pass, history-based deciduous test | 95 | 30 (32%) | 32 | 28 | 30 | 13 |
| Oct 11 pass, seasonal deciduous test | 28 | 14 (50%) | 4 | 9 | 14 | 10 |

Per class with the seasonal test: peak 7 of 18 correct (39%), turning 7 of 10 (70%).

**The seasonal test.** Winter (Jan-Feb 2025 and 2026, five clear passes) and summer (Jun-Aug 2026, six passes) NDVI medians on the same 10 m grid. Deciduous likelihood = (1 − winter/summer) / 0.25, clipped to 0..1, with pixels of summer NDVI ≤ 0.50 set to 0; a pixel needs 0.10 or more to be called turning or peak. By label, the median winter/summer ratio was 1.01 for peak cells, 1.05 turning, 1.11 green deciduous, 1.15 conifer, 1.19 lawn and 1.28 built; the AUC of the ratio for colour cells against conifer, lawn and built was 0.69 to 0.71, against 0.52 to 0.58 for every autumn-only feature in the first analysis. Lawns are greener in January than in a dry Seattle September, which is why a winter-versus-September ratio (tried first) did not work and a winter-versus-summer ratio does.

**Why not a canopy layer.** The City of Seattle tree canopy data was the first choice but stops at the city limit, so it could not be scored here, and its host (services.arcgis.com, data-seattlecitygis.opendata.arcgis.com) is not reachable from the build environment. The 1 m canopy-height map from Meta and WRI (imagery-derived, on AWS Open Data) was tried instead and rejected: over this site it captures only a fraction of the residential crowns the drone image shows (`meta_canopy_vs_drone.png`; correlation with drone vegetation 0.32; cells whose crown fills most of the square have a median canopy fraction of 0.04 in it).

**Cost.** Recall halves. The cells lost are partial crowns sharing a 10 m square with winter-green lawn and shrubs; at 10 m those look evergreen. Five of the 18 full-crown colour cells have a winter/summer ratio above 1 and are missed for the same reason or because of label noise. Most remaining false positives are pixels with no clear winter or summer pass, which fall back to the history.
