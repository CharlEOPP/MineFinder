# 0004: Merge label sources greedily by priority, deduping only across sources

**Decision:** Known-mine labels come from three inventories. We merge them by:
1. **Source priority**: USMIN, then Colorado inactive mines (`co_inactive`), then MRDS. The best source's record sets each label's location and name.
2. **Cross-source dedupe only**: records from *different* sources within **50 m** collapse into one label. Two records from the *same* source are never merged, however close.
3. **Greedy clustering**: process records best-source-first, and each record joins the nearest existing label within 50 m that doesn't already contain its source.

## In plain words

The same old shaft often appears in all three databases, each with a slightly different dot on the map:

- **USMIN** (USGS) was digitised from topographic map symbols, so its dots are usually within a few metres of the real feature.
- **Colorado inactive mines** come from state field inventories and are generally good.
- **MRDS** (USGS, older) often gives a location only to the nearest section or claim, which can be hundreds of metres off.

If we kept all three dots, the model would get three "mines" for one shaft. Two of them would be on empty hillside, which teaches it that empty hillside is a mine. So we collapse them into one and keep the most accurate dot.

**Why never merge within a source?** USMIN might list two shafts 20 m apart because there really are two shafts. The inventory already decided those are separate features, and we shouldn't overrule it.

**Why greedy and not "join everything within 50 m"?** Chaining. Picture shaft A (USMIN), an MRDS dot 40 m east, and shaft B (USMIN) 40 m further east. "Join anything within 50 m" links A–MRDS–B, and suddenly two real shafts become one label. Greedy-by-priority makes A and B their own anchors first, and the MRDS dot attaches to whichever is nearer.

## Formally

- This is **entity resolution** on spatial records, with a **distance threshold** plus a **cannot-link constraint** (same-source pairs). Plain single-linkage clustering (connected components of the "within r" graph) is transitive and violates the constraint. Our greedy anchor assignment isn't transitive, so it respects it.
- Algorithm, with n records sorted by `(priority rank, source_id)`:
  ```
  for r in records:
      candidates = clusters whose anchor is within dedupe_m of r
                   and whose sources don't include r.source
      if candidates: add r.source to the nearest one
      else:          start a new cluster anchored at r
  ```
  Distances are measured to the **anchor**, the best-source point, not to the cluster centroid, so a cluster can't drift toward the coarse MRDS location.
- **Determinism**: sorting by `(rank, source_id)` before clustering makes the output independent of input order. `label_id = "{source}:{source_id}"` of the anchor is stable across runs and shows where each label came from.
- **Complexity**: currently about O(n²), because the STRtree is rebuilt each step. That's fine for one district (hundreds of records) and gets replaced before the statewide Spark run, where Sedona does a spatial join on a grid.

## Why 50 m

It's big enough to catch the normal disagreement between USMIN and the state inventory, which is usually tens of metres. It's small enough that two distinct workings on the same hillside rarely get merged. MRDS records more than 50 m off just stay as their own labels. That's deliberate: a far-off MRDS dot is a *weak* label, and Phase 2's label-quality step will decide whether to keep, snap or drop it. We don't want to hide it in a merge.

## Why this and not the alternatives

| Option | Why not |
|---|---|
| Keep all records | Duplicate positives, and the coarse ones label empty ground |
| Union all within 50 m | Chains distinct features together |
| Average the locations | Pulls the accurate USMIN point toward a coarse MRDS one |
| Large radius (e.g. 200 m) to catch MRDS | Merges real neighbouring workings in a dense district |

## Where it lives

- `src/minefinder/acquire/labels.py` → `SOURCE_PRIORITY`, `normalize`, `merge_labels`
- `config/aoi_ward_goldhill.yaml` → `labels.dedupe_m`
- Tests: `tests/acquire/test_labels.py` covers all of these, including `test_same_source_points_never_merge`, `test_merged_location_comes_from_highest_priority_source` and `test_label_ids_unique_and_deterministic`

## Revisit when

- Real data shows lots of USMIN↔state pairs between 50 and 100 m. Then measure the distribution and pick the threshold from it.
- **Proposed brief change**: record source priority and the cross-source-only rule in `minefinder-brief.md`. This is pending your approval.
