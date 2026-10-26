# 0006: Pick specific label sources, and clean each one's quirks before merging

**Decision:** The three label sources are:

| Slot (`source`) | Dataset | File in `data/raw/` |
|---|---|---|
| `usmin` | USGS Prospect- and Mine-Related Features from Topographic Maps, **ver 10.0**, Colorado files | `usmin/CO-point.shp`, `usmin/CO-poly.shp` |
| `co_inactive` | Colorado Geological Survey **ON-008-04D Abandoned Mine Land Inventory, USFS lands** (field survey, 1991–1999) | `co_aml/…/USFS AMLI.gdb` |
| `mrds` | USGS **Mineral Resources Data System**, national CSV, filtered to Colorado | `mrds/mrds.csv` |

Each one gets its own loader in `acquire/load_labels.py`, which fixes that source's quirks *before* the cross-source merge from [0004](0004-label-merge-and-dedupe.md). All of these quirks were found by measuring the real downloads, not assumed.

## In plain words

Each inventory was made by different people, for different reasons, with different tools. Each one has its own way of recording the same mine twice, or recording things that aren't what we're looking for. If we fed them straight into the merge, the model would learn from those mistakes.

The rule from 0004 is "never merge two records from the *same* source", because an inventory knows its own features. That only holds if the inventory has no internal duplicates. Every one of our three turned out to have them. So each loader cleans up its own source first, and the merge only has to handle disagreements *between* sources.

## What we found, and what each loader does

These numbers are from the Ward / Gold Hill / Caribou area of interest.

### USMIN (573 points + 4 polygons in the AOI)

| Finding | Action | Why |
|---|---|---|
| **No feature ID column** | Derive one: a 12-character SHA-1 of (layer, type, quad, scale, coords) | We need an ID that stays the same between runs. Hashing it means the coordinates can't leak out through an ID in an export. |
| **ver 10.0 never reconciled map scales.** 72 points come from 1:62,500 quads, and **71 of them sit within 41 m of a 1:24,000 point** (median 11 m) | Drop any non-24k symbol within 50 m of a 24k symbol | It's the same shaft drawn on two maps. The 24k map is about 2.6× more detailed, so its symbol is the better location. |
| 3 **gravel pits** in the AOI (and borrow/sand pits statewide) | Drop types containing gravel / borrow / sand pit / clay pit / cinder | Those are road-building aggregate pits, not hard-rock workings. They're exactly the "expected false positive" landforms from the brief, so labelling them as mines would teach the wrong lesson. |

USGS has since published a *Consolidated* release that does the scale merge itself. The loader handles ver 10.0 correctly, so there's no need to re-download. If we ever switch to the consolidated release, the scale filter just finds nothing to drop.

### CGS USFS Abandoned Mine Land Inventory (465 openings, 368 dumps in the AOI)

| Finding | Action | Why |
|---|---|---|
| Already in **EPSG:26913** | Nothing to do | |
| Openings (adits, shafts, prospect holes) and dumps are **separate layers** | Load both as labels, with `feature_type` from `HTYPE` / `PTYPE` | An adit and its dump are different landforms with different shapes in the terrain. Both are useful. |
| **59% of dumps share their opening's exact coordinates** | Drop a dump within 1 m of an opening | The surveyor reused one GPS fix for the whole site. That's one location, not two features. Dumps with their own fix are kept. |
| Feature records have no names, but `Unique_ID` = area ID + 3-digit feature number, and the **area** record has `SITENAME_1` | Join on `Unique_ID[:-3]` | This matched **100%** of openings and dumps statewide. |

### MRDS (722 records in the AOI)

| Finding | Action | Why |
|---|---|---|
| Coordinates are **WGS84** (`mrds.met`: "Horizontal_Datum_Name: World Geodetic System 1984") | Read as EPSG:4326 | If we'd assumed NAD27, a common default for old USGS data, every point would shift by tens of metres. |
| **198 records sit on an identical coordinate** with another record | Keep the lowest `dep_id` at each coordinate | Several mines geocoded to the same section or claim centroid. That's one coarse location, and none of the stacked copies is actually on its mine. |
| `feature_type` comes from `dev_stat` (past producer, occurrence, …) | Keep it, but note it's *development status*, not a landform type | It's still useful later for weighting labels (a past producer has workings, an occurrence may not). |

## What the cross-source distances told us

Distance from each record to the nearest USMIN 1:24k point:

| From | median | ≤ 50 m | ≤ 250 m |
|---|---|---|---|
| CGS openings | 850 m | 14% | 21% |
| CGS dumps | 887 m | 15% | 21% |
| MRDS | 310 m | 19% | 40% |

Two conclusions:

1. **There's no datum shift between CGS and USMIN.** For the 86 pairs within 150 m, the median offset is (−10 m, 0 m), with about 35 m scatter in both directions. A NAD27/NAD83 mix-up would show up as a consistent shift of tens of metres in one direction. Instead it's random scatter, as you'd expect from hand-held GPS and map symbols.
2. **CGS and USMIN are mostly complementary, not overlapping.** Both cover most of the AOI (36 of the 2 km cells contain both), but they rarely record the *same* features. A likely reason is that CGS only surveyed **Forest Service land**, while much of the district's core is **patented private claims**, which only the topo maps cover. This is a hypothesis for now. The Phase 5 land-status layer will confirm or reject it. Either way, adding CGS roughly **doubles** the field-verified labels.

So the 50 m threshold from 0004 holds up. Real duplicates cluster well inside it: the USMIN scale duplicates have a median of 11 m. Records further apart than 50 m are mostly genuinely different features. MRDS is the exception: its records are often far from the real feature *and* near other ones. That's why it gets [weak-label handling](#revisit-when) rather than a bigger radius.

## Result (current AOI)

`python -m minefinder.acquire.load_labels` writes `data/interim/labels.gpkg`, with **1,370 labels**:

| sources | count |
|---|---|
| co_inactive only | 448 |
| mrds only | 409 |
| usmin only | 345 |
| mrds + usmin | 91 |
| co_inactive + usmin | 53 |
| co_inactive + mrds | 12 |
| all three | 12 |

By split: **train 1,090 · test 274 · buffer 6** (buffer labels are discarded, per [0003](0003-spatial-holdout-and-buffer.md)). The `split` column is set at build time, so no later step can forget to check it.

## Formally

- This is **per-source data cleaning ahead of entity resolution**. Each loader enforces the merge's precondition, "no within-source duplicates", using the source's own structure (map scale, a shared GPS fix, an identical geocode), not a generic distance rule.
- The USMIN ID is a **content hash**: deterministic, collision-checked (exact duplicate records get a `-1`, `-2` suffix), and **non-reversible**, so it can't leak coordinates.
- The thresholds are 50 m for scale duplicates (the same as `dedupe_m`, and above the largest observed pair at 41 m) and 1 m for a shared GPS fix (well above float noise, well below any real feature spacing).

## Why this and not the alternatives

| Option | Why not |
|---|---|
| Download the USMIN Consolidated release instead | Would work, but ver 10.0 is already downloaded and the scale filter makes the two equivalent for us |
| USGS OFR 94-579 as the Colorado source | A desk compilation from 1994. CGS is a field survey, which is better location quality for the same slot. Keep OFR 94-579 as a fallback |
| Treat USMIN 1:62.5k as a separate lower-priority "source" | Same result, but it adds a fake source to `SOURCE_PRIORITY`. A loader-level filter is simpler |
| Drop MRDS entirely | 409 MRDS-only labels may point at real, otherwise unlabelled workings, so they're worth keeping as *weak* labels |
| Keep gravel pits as positives | They're the false-positive class from the brief |

## Where it lives

- `src/minefinder/acquire/load_labels.py`: `load_usmin`, `load_co_aml`, `load_mrds`, `build_labels`, the CLI
- `src/minefinder/acquire/labels.py`: `normalize(..., type_col=)` and `feature_type` in the merged schema
- Tests: `tests/acquire/test_load_labels.py`, which has one test per finding above, using synthetic files in each real format, plus a real-data smoke test that skips when `data/raw/` is missing

## Revisit when

- **Phase 2 label quality.** MRDS-only labels (409, about 30%) are weak. Decide whether to snap them to the nearest LRM landform within about 250 m, down-weight them, or drop them. Measure on the train region only.
- **The land-status layer arrives (Phase 5).** Test the "CGS = USFS land, USMIN = patented core" hypothesis.
- **Statewide.** The CGS inventory only covers USFS land. BLM land needs its own AML source (the BLM AML inventory).
