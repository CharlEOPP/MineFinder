# 0010: A rule-based baseline, and how every detector gets scored

**Decision:** Phase 3 delivers two things:

1. **An evaluation protocol** (`detect/evaluate.py`) that every detector from now on is scored with: **precision@N with a 40 m match radius and one-to-one matching**, always reported next to the **random-point precision** and **lift** (precision ÷ random). It's computed against all labels and against USMIN + CGS only.
2. **A rule-based baseline** (`detect/baseline.py`): pits from `lrm_small`, dumps from `lrm_large`, rejected by shape (roads, ditches, ridges) and by the seam mask, then ranked by prominence. Thresholds are tuned **on the train region only**, and the held-out test block is scored **only with `--final`, which is logged every time** in `docs/test_evaluations.csv`.

**Label snapping is deferred to Phase 4, and test labels are never snapped.** (This changes the plan in 0008 and 0009.)

## In plain words

**How we score a detector.** The detector hands us a ranked list of "possible mines". We take the top 50 and ask how many are within 40 m of a known mine. That fraction is **precision@50**. It's the number that matters for field work: "if I check the top 50, how many are real?"

Three details keep it honest:
- **One find per mine.** If the detector puts five candidates on the same dump, that's one hit and four misses. Otherwise a detector could score well by repeating itself.
- **Compare against luck.** Mines are dense in this district, so even a random dot is sometimes within 40 m of one. We measure that chance by scattering 20,000 random points, and report **lift** = precision ÷ chance. A precision of 0.12 means little alone, but "16× better than random" is meaningful.
- **It's a lower bound.** A candidate on a real but *unrecorded* mine counts as a miss, because it isn't in any inventory. Those are exactly the finds we're after. So real precision is at least this number, and Phase 5's manual review is what reveals the rest.

**Why 40 m:** the CGS labels sit tens of metres from their features (0006, 0008), so a tighter radius would mark real finds as misses because the *label* is off. A looser one lets neighbouring workings "help" each other. We also report 25 m and 60 m to show how sensitive the result is.

**Why no snapping for evaluation:** snapping moves each label to the nearest pit or dump, which is precisely what the detector finds. Snapping test labels would pull them toward the detector's answers and inflate its score, which amounts to leaking the detector into the answer key. The match radius handles label error instead. Snapping is still useful, but for building *training* targets in Phase 4, on train labels only.

**How the baseline works.** Mining leaves two kinds of bumps:
- **Pits** (shaft collars, prospect pits): a low spot in the 21 m local relief. We find each connected low patch and place the candidate at its deepest pixel.
- **Dumps** (waste rock, tailings): a high spot in the 61 m local relief, placed at its highest pixel.

Then we throw out look-alikes **by shape**:
- **Roads, ditches, ridges, gullies** are long. Anything wider than 120 m, or more than about 2–4× longer than it is wide, is rejected.
- **Tiny blobs** (under 4 m²) are noise.
- **Anything in the seam mask** is unreliable terrain (0009).
- **"Dumps" taller than about 4 m** are natural knobs. On the train cell, mounds taller than 4 m had a 0% hit rate.

What's left is ranked by **prominence**: pit depth, or dump height × a weight so the two kinds share one list.

**Tuning without cheating.** The thresholds above are chosen by trying 486 combinations and keeping the one with the best precision@50 **on the train region**. Test candidates and test labels are removed *before* tuning starts, and a test proves that changing the test labels can't change the chosen parameters. The test block is scored once, at the end, with `--final`, and every such run appends a line (time, parameters, scores) to `docs/test_evaluations.csv`, which is committed. If we ever peek repeatedly, git history will show it.

## Formally

- **Matching:** candidates sorted by score descending. Each takes its nearest unclaimed label within r (an STRtree `dwithin` query). precision@N = hits in the top N ÷ **N**, not ÷ the number returned, so a detector can't look better by returning fewer than N candidates.
- **Random precision:** the fraction of 20,000 uniform points in the evaluation region within r of any label. It's the expected precision of a detector with no skill. Lift = P@N ÷ random.
- **Components:** 8-connected via `scipy.ndimage.label` on `lrm_small < −0.2 m` (pits) and `lrm_large > +0.3 m` (dumps), the loose floors. Per component: area, bbox extent, the extreme-value location, and **elongation** = √(λ_max/λ_min) of the pixel-coordinate covariance (1 = round).
- **Blockwise:** 2048 px blocks with a 128 px halo. A component is kept only by the block whose core contains its location, and components wider than `max_extent_m` (120) are dropped everywhere, so results equal a whole-raster run exactly (tested).
- **Extract once, tune many:** extraction stores every loose candidate with its measurements (`outputs/baseline_candidates.gpkg`, which is **gitignored because it holds exact coordinates**). Tuning and scoring only filter and sort that table.
- **Tie-break:** with N = 50, many settings tie. Ties go to the best p@100, then to grid order.
- **Halo guard** (optional, tuned): drop a pit that has a strong dump within 12 m but isn't on it, since that's likely the 21 m window's halo (0009). On the train cell, tuning chose **off**, because the depth threshold already removes halo artifacts.

## Dev result (one train-region cell, x45y444, 2026-10-26)

Tuned parameters: `pit_min_depth 0.3, dump_min_height 0.5, dump_max_height 4.0, max_elongation 2.0, dump_weight 0.25, halo_guard off`.

| Labels | P@25 | P@50 | P@100 | Random (40 m) | Lift @50 |
|---|---|---|---|---|---|
| all (161) | 0.08 | **0.12** | 0.08 | 0.0075 | **16×** |
| USMIN + CGS (105) | 0.08 | **0.12** | 0.07 | 0.0051 | **24×** |

**What the data said while building it:**
- **Deeper pits are more likely to be known mines:** a 0.15% hit rate at 0.3–0.5 m deep, 2.8% at 1–2 m, 7.9% at 2–3 m, and 15% beyond 3 m. Ranking pits by depth is right.
- **Dumps peak at 2–4 m high** (1.2% hit rate). Taller "dumps" scored 0%, which is why `dump_max_height` was added.
- **A dump within 15 m doubles a pit's hit rate** (0.61% vs 0.29%), but only 5% of pits have one. That's a candidate feature for Phase 4, not a baseline rule.
- **The grid is flat at the top** (36 of 486 combinations tie at 0.12). The baseline has hit the ceiling of "rank by one number". The top 50 are all deep pits, and 44 of them aren't known mines. Some are rock hollows; some may be unrecorded workings, which is the point of Phase 5.

**This dev number is not the result.** It's one cell, entirely in the train region. The real baseline number is the **full-AOI `--final` run on the held-out Caribou block**, recorded in `docs/test_evaluations.csv`. That's the number Phase 4 has to beat.

## Official baseline (full AOI, 2026-10-26)

Tuned on the train region (`config/baseline_params.yaml`): `pit_min_depth 0.3, dump_min_height 2.0, dump_max_height 4.0, max_elongation 2.0, dump_weight 1.0, halo_guard on`. Only 3 of 486 settings tied at the top.

| Region | Labels | P@25 | **P@50** | P@100 | Random (40 m) | Lift @50 |
|---|---|---|---|---|---|---|
| Train (dev) | all (1,090) | 0.04 | 0.12 | 0.07 | 0.012 | 10.0× |
| Train (dev) | USMIN+CGS (726) | 0.04 | 0.08 | 0.05 | 0.008 | 10.3× |
| **Test (`--final`, once)** | **all (274)** | **0.16** | **0.08** | **0.07** | **0.022** | **3.7×** |
| **Test (`--final`, once)** | USMIN+CGS (231) | 0.16 | 0.08 | 0.07 | 0.018 | 4.5× |

**The number Phase 4 must beat: precision@50 = 0.08 at 40 m (4 of 50) on the held-out Caribou block, 3.7× random.** It's logged in `docs/test_evaluations.csv`.

**How to read it:**
- **It's small-sample.** 4 hits out of 50 has a 95% interval of roughly **0.03–0.19** (Wilson). The train value, 6/50, is 0.06–0.24. The intervals overlap, so "test is worse than train" is plausible but not proven. Phase 4 comparisons should look at P@25, P@50 and P@100 together, not one cell.
- **The lift drop (10× → 3.7×) is mostly the denominator.** Caribou is denser in known mines (random precision 2.2% vs 1.2%), so chance does better there. Precision itself fell from 0.12 to 0.08.
- **The ranking front-loads its best guesses.** P@25 = 0.16 (4 of 25) on test: all 4 hits are in the top 25. Beyond that, the next 75 add only 3 more.
- **It's a lower bound.** The 46 non-hits are Phase 5 review material: some will be rock hollows, and some may be unrecorded workings.

## Bug caught on Charlie's PC

`test_blockwise_extraction_equals_whole` passed here but failed on Windows: one dump landed 1 m apart between the whole-raster and blockwise runs. **Cause:** LRM is stored in whole centimetres, so a dump's top is often a *plateau* of tied pixels, and `scipy.ndimage.maximum_position` / `minimum_position` don't promise which tied pixel they return. A new test showed that even the scipy here doesn't return the first one, so the earlier agreement was luck. **Fix:** find the extreme value per component, then take the **first tied pixel in row-major order** (top-most, then left-most). That's explicit, version-independent and translation-invariant, so every block agrees. The dev result was unchanged (P@50 = 0.12). This is exactly what a blockwise = whole test is for: it turned a silent, platform-dependent off-by-one into a red test.

## Why this and not the alternatives

| Option | Why not |
|---|---|
| Precision over *all* candidates, or F1 | Field effort is spent on the top of the list, so precision@N is the decision-relevant number (it's in the brief) |
| Many-to-one matching | Lets duplicates inflate precision |
| Snapped test labels | Leaks detector-like features into the answer key |
| Tuning on the test block | Turns the test block into a train block, so the headline number would be fiction |
| Learned weights in the baseline | That's Phase 4. The baseline should be transparent rules anyone can audit |

## Where it lives

- `src/minefinder/detect/evaluate.py`: `match_candidates`, `precision_at_n`, `random_precision`, `evaluate`
- `src/minefinder/detect/baseline.py`: `extract_candidates`, `score_candidates`, `tune`, `evaluate_split`, and a CLI with `extract | tune | evaluate [--final]`
- `config/baseline_params.yaml` (written by `tune`; committed), `docs/test_evaluations.csv` (written by `--final`; committed)
- Tests: `tests/detect/test_evaluate.py` and `tests/detect/test_baseline.py`. They cover synthetic pits, dumps, a road fill and a ditch; blockwise = whole; tuning blind to test labels; test scoring that needs `--final` and is logged; the tie-break; and the height cap.

## Revisit when

- **Phase 4:** a learned model gets the same `evaluate()` and the same `--final` log. Features worth learning: dump-within-15 m, depth × area interplay, slope context (adits sit on slopes with a bench below).
- **Phase 5:** the 44 "misses" in a top-50 are review material. Every candidate must pass the land-status filter and the MLRS active-claim check before any visit, and no exact coordinates go in anything published.
- **Label quality improves** (snapping in Phase 4, or better CGS locations): re-run `--final` with the same parameters, so the comparison stays apples to apples.
