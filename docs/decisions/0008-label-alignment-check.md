# 0008: Measure label alignment against random points, with a shift search only where it means something

**Decision:** `review/alignment.py` measures whether known-mine labels sit on their landforms, group by group (anchor source × landform class). It compares each group against **random points in the same DEM**. A **shift search** (±40 m) estimates any systematic offset, and it's significance-tested against the same search run on random points. The search runs only for classes where "nearest strong low/high" really means "the feature": **pit-like** (looks for low ground) and **dump** (looks for high ground). By default it uses **train-split labels only**.

## In plain words

"Does this dot sit on a mine?" is hard to answer by eye at map scale. So we ask the terrain instead.

1. **Local relief** = elevation minus the average elevation of the surrounding 31 × 31 m. A hillside becomes ~0, a shaft collar or prospect pit becomes clearly negative, and a waste dump becomes positive. (This is the first piece of the Phase 2 "local relief model".)
2. **Score each label:** how deep is the deepest low within 4 m of it? The 4 m tolerance lets a dot on the pit's *rim* still count.
3. **Compare against chance.** Random points also land near lows sometimes, because terrain is bumpy. So we score 3,000 random points and ask how often a random set of the same size does as well as the labels. That fraction is the **p-value**. A small p (< 0.05) means "better than luck".
4. **Look for a systematic offset.** Slide all labels in a group together by every shift up to ±40 m and find where they line up best. A CRS or datum mistake would show up here as one consistent shift.

**The catch in step 4:** if you try 6,561 shifts, one of them always looks best, even for random dots. So the best shift is tested the same way: we run the identical search on random point sets and ask how often *they* reach as good a "best". That's `p_best_shift`. Without it, every group would appear to have an offset.

## Formally

- **Statistic:** the group mean of `pit_depth(sign × LRM, r=4)` at the label (`score_at_label_m`). sign = −1 for dumps, so mounds are scored as height.
- **Null:** a bootstrap of group-size means from the random-point pool (2,000 draws), giving a one-sided p.
- **Shift estimate:** the argmax over (dx, dy) ∈ [−40, 40]² of the group-mean **exact-pixel** relief (−sign × LRM). It deliberately does *not* use the 4 m tolerant depth. That score is flat-topped within 4 m of a feature, so ties resolve to the first index scanned, which biased early versions by (−2, −2) m. The tests caught that.
- **Shift significance:** the null is the max over the same 81 × 81 search for random groups of size n (200 draws). This is a max-statistic permutation test, which corrects for searching many shifts at once (the "look-elsewhere" effect).
- **Windows:** each point reads a 123 × 123 px DEM window. Points whose central region has any nodata are **skipped** and counted (`n_skipped`), because a box mean on a slope with missing pixels is biased (see `test_local_relief_ignores_nan_neighbours`). This keeps memory tiny regardless of mosaic size.
- **No shift search for adits or "other":** an adit is a notch on a slope, so the strongest nearby low is usually the **gully** below it. On the real tile, adits "aligned" 29 m away with p = 0.000, and that turned out to be the metric finding streams. MRDS "other" mixes development-status records with no consistent landform.
- **Train split by default (`--split train`):** we'll tune label snapping from these numbers in Phase 2. Tuning on the test block would be leakage ([0003](0003-spatial-holdout-and-buffer.md)).
- **No coordinates in the output:** the report is aggregate statistics only, so it's safe to share ([0005](0005-repo-layout-and-what-never-gets-committed.md)). A test enforces this.

## Validation (synthetic DEM with 60 known pits / mounds)

| Case | Expected | Result |
|---|---|---|
| Labels exactly on pits | far above random, shift (0, 0) | ✔ score 2.4 m vs 0.14 m random, shift (0, 0) |
| Labels 12 m E, 8 m S of pits | shift (−12, +8), significant | ✔ within 1 m, p < 0.01 |
| Labels 10 m W, 6 m N of mounds (dumps) | shift (+10, −6), significant | ✔ |
| Random labels | not significant, at label *or* best shift | ✔ |

## First real result (one cell, x45y444, all splits)

| Group | n | Score at label | Random | p at label | Best shift | p best shift |
|---|---|---|---|---|---|---|
| CGS pit-like | 45 | 0.19 m | 0.31 | 0.99 | (+25, −4) | 0.64 |
| CGS adit | 44 | 0.47 m | 0.31 | 0.015 | — | — |
| CGS dump | 9 | 0.32 m | 0.33 | 0.45 | (+16, −17) | 0.42 |
| MRDS other | 55 | 0.42 m | 0.31 | 0.04 | — | — |

**No group has a significant offset**, which confirms the Phase 1 "no CRS or datum error" result. CGS shaft/pit labels don't sit on pits any more often than random points do. That's consistent with the GPS-era location error discussed in [0006](0006-label-sources-and-loader-cleaning.md). One tile has too few USMIN labels to report. The full-mosaic run will fix that.

## Full-mosaic result (train split, 2026-10-26)

| Group | n | Score at label | Random | p at label | Best shift | p best shift |
|---|---|---|---|---|---|---|
| USMIN pit-like | 155 | 0.54 m | 0.40 | 0.007 | (−34, +31) | 0.91 |
| USMIN adit | 106 | 1.20 m | 0.40 | < 0.001 | — | — |
| USMIN other | 20 | 0.97 m | 0.40 | 0.001 | — | — |
| CGS pit-like | 157 | 0.22 m | 0.40 | 1.00 | (−37, −40) | 0.96 |
| CGS adit | 212 | 0.58 m | 0.40 | < 0.001 | — | — |
| CGS dump | 45 | 0.38 m | 0.45 | 0.68 | (−36, −10) | 0.75 |
| MRDS other | 357 | 0.59 m | 0.40 | < 0.001 | — | — |

- **There's no systematic offset in any group** (every p_best_shift ≥ 0.75). CRS and datum are confirmed on about 1,100 labels.
- **USMIN adits sit on clear notches**, at 3× the random relief. USMIN pits are better than random but often a few metres off a pit.
- **CGS pits are worse than random** (p = 1.00). Either the 1990s GPS error is large compared to pit size, or many have been backfilled, or both. CGS adits are fine.
- **Decision: snap labels later, not now.** Snapping should pick "the pit" or "the notch", which needs the Phase 2 derivatives and a feature-specific score. It becomes part of building training targets in Phase 3, with its radius chosen from these numbers on the train split.

## Where it lives

- `src/minefinder/terrain/derivatives.py`: `local_relief`, `pit_depth` (numpy only, with integral-image box means)
- `src/minefinder/review/alignment.py`: `label_groups`, `alignment_report`, and the CLI (writes `data/processed/alignment.csv`)
- Tests: `tests/terrain/test_derivatives.py` and `tests/review/test_alignment.py`

## Revisit when

- **Adits:** score them with a *notch + bench* signature (a slope break with a flat bench downhill) instead of "nearest low". That's the Phase 2 detector's job, and the alignment check can reuse it.
- **Label snapping (Phase 2):** if USMIN shows a real but small offset on the full mosaic, snapping radius and direction come from this report (train split only).
- **Survey seams:** random points that land on a seam line score as weak lows. That's negligible now. Recheck after seam handling.
