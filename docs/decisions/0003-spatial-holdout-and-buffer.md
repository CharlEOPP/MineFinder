# 0003: Use a geographic test block plus a buffer, and discard anything that straddles them

**Decision:** The test set is a fixed **geographic block**, the Caribou corner of the AOI, not a random sample. A **256 m buffer** around it belongs to neither train nor test. Any label or image chip that isn't *entirely* inside the test block or *entirely* inside the train region is labelled `buffer` and thrown away.

## In plain words

Imagine you split your data by flipping a coin for every 256 m image chip. Two neighbouring chips overlap the same waste-rock dump, and one goes to train and the other to test. The model memorises "this exact lumpy texture = mine" from the train chip, then "detects" it in the test chip. Your test score looks great, but the model was just recognising a place it had already seen. That's **leakage**.

Mines cluster too. The Ward district has dozens of workings along the same veins, so a random split puts siblings of every test mine into training.

The fix is to wall off a whole region. Everything in the Caribou block is test, and the model never sees it during training *or tuning*. That's the honest version of the question we actually care about: "Will this find mines in places it hasn't seen?"

The buffer handles the edge. A chip just outside the wall can still see across it, so we leave a gap of no-man's-land at least as wide as a chip plus its context window.

## Formally

- **Spatial autocorrelation**: nearby observations are more similar than distant ones (Tobler's first law). Random k-fold CV assumes independent and identically distributed samples, so on spatial data it gives **optimistically biased** error estimates.
- **Spatial blocking / block cross-validation** (Roberts et al., 2017, *Ecography*): hold out contiguous regions. Add a **buffer** wider than the autocorrelation range, or the receptive field of your model, whichever is larger.
- `split_of(geom)` assigns a label in a fixed order:
  1. not within the extent → `outside`
  2. within the holdout → `test`
  3. within (extent − holdout.buffer(b)) → `train`
  4. otherwise → `buffer` (discarded)

  Using `within`, not `intersects`, is what makes a straddling chip fall through to `buffer`.
- The buffer size `b ≥ chip_size + context` (256 m now). If Phase 3/4 uses larger chips, `b` must grow to match.

## Why this and not the alternatives

| Option | Why not |
|---|---|
| Random train/test split | Leaks through overlapping chips and clustered mines, so scores look much better than they really are |
| Holdout with no buffer | Chips on the boundary still see across it |
| Keep straddling chips in train | They contain test pixels, which is leakage by definition |
| Several spatial folds now | Worth doing later (Phase 4 eval); one block is enough to build the pipeline honestly |

**Why Caribou:** it's a corner of the AOI, so removing it doesn't cut the training region in two. It also has a different mix of workings from Ward, so it tests whether the model *generalises* rather than only interpolates.

## The rule this protects

From the project instructions: **no result claims without precision@N on the held-out spatial block.** This decision is what makes that number mean something. Tuning thresholds while looking at the Caribou block would spend it. If we need a tuning set, we carve a *second* block out of the train region.

## Where it lives

- `config/aoi_ward_goldhill.yaml` → `holdout`, `holdout_buffer_m`
- `src/minefinder/aoi.py` → `AOI.train_region`, `AOI.split_of`
- Tests: `test_train_region_excludes_holdout_and_buffer`, `test_train_region_area`, `test_split_of` (the straddling-chip cases)

## Revisit when

- Chip size or context window changes, since the buffer has to grow with it.
- We add a validation block for tuning (Phase 4).
- **Proposed brief change**: record the 256 m buffer and the discard-on-straddle rule in `minefinder-brief.md`. This is pending your approval.
