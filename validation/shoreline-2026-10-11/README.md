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
