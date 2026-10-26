# MineFinder docs

This folder records **why** the code is the way it is. Code says what, tests say what must stay true, these docs say why.

## How to read them

- `decisions/` holds one file per decision, numbered in the order they were made. Each one follows the same shape:
  1. **Decision**: the one-line answer.
  2. **In plain words**: the idea with a concrete example, no jargon.
  3. **Formally**: the precise version and the terms you'd search for.
  4. **Why this and not the alternatives**: what else was on the table and why it lost.
  5. **Where it lives**: the files and tests that enforce it.
  6. **Revisit when**: the condition that would make us change it.
- `sessions/` holds a short log per working session: what got done, which decisions were made, the next concrete step.

A decision is never silently changed. If we change one, the old file gets a "Superseded by NNNN" line at the top and a new numbered file explains the change.

## Decision index

| # | Decision | Phase |
|---|----------|-------|
| [0001](decisions/0001-environment-conda-forge.md) | Use a conda-forge environment, installed with Miniforge | 0 |
| [0002](decisions/0002-crs-epsg-26913.md) | Use EPSG:26913 everywhere and reject anything else at the boundary | 0 |
| [0003](decisions/0003-spatial-holdout-and-buffer.md) | Use a geographic test block plus a buffer, and discard anything that straddles them | 0–1 |
| [0004](decisions/0004-label-merge-and-dedupe.md) | Merge label sources greedily by priority, deduping only across sources | 1 |
| [0005](decisions/0005-repo-layout-and-what-never-gets-committed.md) | Use a src-layout package, keep data out of git, and keep candidate coordinates private | 0 |
| [0006](decisions/0006-label-sources-and-loader-cleaning.md) | Pick specific label sources, and clean each one's quirks before merging | 1 |
| [0007](decisions/0007-dem-acquisition-and-mosaic.md) | Get the DEM from the TNM API, and mosaic newest survey first without ever resampling | 1 |
| [0008](decisions/0008-label-alignment-check.md) | Measure label alignment against random points, with a shift search only where it means something | 2 |

## Phase tracker

| Phase | Status |
|-------|--------|
| 0: setup | ✅ done: env, repo, AOI config, test suite |
| 1: labels + DEM | ✅ done: 1,370 labels, 1 m DEM mosaic, no CRS or datum offset. CGS labels sit tens of metres from their features (see 0006) |
| 2: terrain | 🟡 alignment tool built and tested (0008); next: run it on the full mosaic, then derivatives, LRM and seam handling |
| 3–6 | not started |
