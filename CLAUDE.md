# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

A **learning project**, not a production codebase. The owner is a second-year
Molecular and Cell Biology major building pandas/numpy fluency and a portfolio
piece. The analysis question — does tumor mutational burden predict survival on
immune checkpoint inhibitors — is the vehicle, not the point.

Two consequences that override normal defaults:

1. **The human writes the substantive code.** See the next section. This is the
   single most important thing in this file.
2. **Portfolio quality is a requirement.** Reproducibility, tests, and an honest
   limitations section are deliverables, not polish.

## Do not write the code marked `YOU WRITE THIS`

The implementation plan marks certain functions `YOU WRITE THIS`. Do not
implement them — not when asked to "just fix it", not when the human is stuck,
not when writing it would be faster. Code Claude wrote is worthless for both of
this project's goals: it cannot be claimed in an interview, and nothing is
learned from it.

When the human is stuck, **narrow the gap instead of closing it**: ask what they
tried, name the relevant method, explain the concept, or give a small analogous
example over different data.

Exception: if they explicitly ask for an implementation after attempting it, write
it — and say plainly that they should rewrite it themselves afterward.

**Fine to write outright:** scaffolding, config, `__init__.py`, Makefile,
requirements, and the test suites the plan supplies verbatim. The constraint
covers code embodying a real decision or a transferable pandas idiom.

## Commands

```bash
make data        # fetch from cBioPortal into data/raw/ (no-op once cached)
make process     # build data/processed/analysis_table.csv
make figures     # regenerate figures/
make test        # pytest -v
make all         # data -> process -> figures
make clean       # drop derived data and figures (never touches data/raw/)
```

Run a single test file, or a single test:

```bash
pytest tests/test_fetch.py -v
pytest tests/test_tmb.py::test_percentile_is_computed_within_group_not_globally -v
```

Always run `pytest` from the repository root — `src` is only importable from
there. Tests must never require network access; the fetch tests monkeypatch
`requests.get` and `time.sleep`.

Python is Anaconda's `/opt/anaconda3/bin/python` (3.13.9). pandas, numpy,
matplotlib, seaborn, scipy, scikit-learn, and pytest are already installed.
`lifelines` is **not** installed and is not needed until Phase 2.

## Architecture

Data flows one direction. No module imports anything downstream of it:

```
fetch -> clean -> features -> {analysis, models} -> viz
```

| Module | Responsibility | Deliberately does not know |
| --- | --- | --- |
| `src/data/fetch.py` | HTTP, retry, on-disk caching | what any field means |
| `src/data/clean.py` | pivot long→wide, dtype coercion, join | TMB, survival semantics |
| `src/features/tmb.py` | pure DataFrame→DataFrame transforms | plotting, modeling |
| `src/analysis/survival.py` | KM, log-rank, Cox (Phase 2) | figure styling |
| `src/models/train.py` | estimators and evaluation (Phase 3) | figure styling |
| `src/viz/plots.py` | figures only | any analysis logic |

**Notebooks contain narrative and figures, never transformation logic.** Anything
worth testing belongs in `src/` and gets imported. This is what makes the repo
reviewable and the transformations testable; preserve it.

Every module under `src/` has a matching `tests/test_<module>.py`.

## Verified facts — assert, do not re-derive

These were checked against the live API on 2026-07-22 and are asserted in
`clean.py`. If a run contradicts one, the pipeline should stop rather than
silently analyze the wrong cohort.

| Quantity | Value |
| --- | --- |
| PATIENT endpoint records | 9,966 |
| SAMPLE endpoint records | 24,076 |
| Patients / samples after pivot | 1,661 / 1,661 |
| Cancer types | 11 |
| Deceased / living | 832 / 829 |
| TMB min / median / p90 / max | 0.00 / 5.87 / 25.45 / 207.5 |
| Gene panels | IMPACT341 (230), IMPACT410 (1,001), IMPACT468 (430) |

## Data sources

Primary: cBioPortal REST API, study `tmb_mskcc_2018` (Samstein et al.,
*Nature Genetics* 2019). No API key.

```
https://www.cbioportal.org/api/studies/tmb_mskcc_2018/clinical-data?clinicalDataType=PATIENT&pageSize=200000
```

**Dead ends already investigated — do not retry:**

- `olink.com/mgh-covid-study/` — HTTP 404, permanently gone; the Wayback capture
  contains no data-file links
- `cbioportal-datahub.s3.amazonaws.com/tmb_mskcc_2018.tar.gz` — HTTP 403; use the
  REST API
- Liu melanoma cohort is dbGaP controlled-access; Gide is EGA raw reads

Phase 4 extension only: GEO `GSE91061` (Riaz 2017), 109 RNA-seq samples from 65
patients. Too small for supervised learning — 51 pre-treatment samples against
~20,000 genes — which is why it is not the spine.

## Domain gotchas

- **Panel heterogeneity is a real confound.** Three panels of differing size
  (341/410/468 genes) make raw TMB non-comparable. The correction is ranking TMB
  *within cancer type* (`groupby().rank(pct=True)`), which is also the project's
  central pandas lesson. Never compare raw TMB across panels.
- **Survival data is censored.** 829 patients are alive *at last contact*, not
  cured. Treating censored patients as survivors biases every estimate.
- **Everything from the API is a string,** including `OS_MONTHS` and
  `TMB_NONSYNONYMOUS`. `OS_STATUS` is the composite `"1:DECEASED"`.
- **TMB is severely right-skewed** and contains exact zeros, so `log1p` rather
  than `log`.
- **The expected result is a weak model.** TMB genuinely explains only part of
  immunotherapy response. Do not tune toward an impressive number; reporting the
  modest effect honestly is the objective.

## Planning workflow

`docs/superpowers/specs/` holds the approved design; `docs/superpowers/plans/`
holds implementation plans.

Plans are written **one phase at a time**, immediately before that phase starts —
never all at once. What the human finds difficult in one phase should reshape the
scope of the next. Phases: 0 setup → 1 wrangling/EDA → 2 survival statistics →
3 ML → 4 optional Riaz extension.

Commit after every task, using the message the plan specifies. The commit history
is itself a deliverable and should show incremental work.

## Writing shell commands for this user

Do not use backslash line-continuations in commands intended to be pasted into a
terminal — they get mangled and the fragments run as separate commands, which has
already caused one broken setup here. Emit self-contained one-liners instead.

This is macOS: `cat -t` (not GNU `cat -A`) to reveal tab characters, which matters
because the Makefile requires real tabs.
