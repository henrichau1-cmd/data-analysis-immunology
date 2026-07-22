# Design: Tumor Mutational Burden and Immunotherapy Survival

**Date:** 2026-07-22
**Status:** Approved
**Repository name:** `immunotherapy-tmb-survival`

## Purpose

A structured, portfolio-quality data analysis project that teaches pandas, numpy,
scipy, and scikit-learn through a real cancer immunotherapy question.

The author is a second-year Molecular and Cell Biology major with a Cognitive
Science minor, planning to take a machine learning course. They can write Python
scripts and have briefly used pandas, but have never completed an end-to-end
analysis. The project must therefore teach data-analysis idioms rather than
Python fundamentals.

Two goals, both binding:

1. **Learning.** Build genuine fluency with the standard Python data stack.
2. **Portfolio.** Produce a public GitHub repository suitable for a resume.

The second goal constrains the first: code quality, reproducibility, and honest
reporting are requirements, not optional extras.

## Research question

**Does tumor mutational burden predict survival in patients treated with immune
checkpoint inhibitors, and does the relationship hold across cancer types?**

The biological rationale is mechanistic and within the author's coursework: more
somatic mutations produce more neoantigens, which give T cells more to recognize,
which should make checkpoint blockade more effective.

This replicates the central finding of Samstein et al., *Nature Genetics* 2019.
Replication is deliberate — a known-correct answer exists, so the author can tell
whether their pipeline is right.

## Data

### Primary dataset (the spine)

**MSK-IMPACT immunotherapy cohort** — Samstein et al., *Nature Genetics* 2019.
Accessed via the cBioPortal REST API, study ID `tmb_mskcc_2018`.

Verified 2026-07-22 by direct API query:

| Property | Verified value |
| --- | --- |
| Patients | 1,661 |
| Samples | 1,661 |
| Cancer types | 11 |
| Survival status | 832 deceased / 829 living |
| Median overall survival | 11.0 months |
| Drug classes | PD-1/PD-L1 (1,307), Combo (255), CTLA-4 (99) |
| Sex | 1,034 male / 627 female |
| TMB (nonsynonymous) | min 0.00, median 5.87, p90 25.45, max 207.5 |
| Gene panels | IMPACT341 (230), IMPACT410 (1,001), IMPACT468 (430) |
| Molecular profiles | mutations, structural variants |

Largest cancer types: NSCLC (350), Melanoma (320), Bladder (215), Renal Cell
(151), Head and Neck (139), Esophagogastric (126), Glioma (117), Colorectal (110),
Cancer of Unknown Primary (88), Breast (44).

Access endpoints (no API key required):

- `https://www.cbioportal.org/api/studies/tmb_mskcc_2018/clinical-data?clinicalDataType=PATIENT&pageSize=100000`
- `https://www.cbioportal.org/api/studies/tmb_mskcc_2018/clinical-data?clinicalDataType=SAMPLE&pageSize=200000`

The bulk download at `cbioportal-datahub.s3.amazonaws.com/tmb_mskcc_2018.tar.gz`
returns HTTP 403 and must not be used.

### Extension dataset (Phase 4, optional)

**Riaz et al. 2017 melanoma cohort** — GEO accession `GSE91061`.

Verified 2026-07-22 via GEO FTP: 109 RNA-seq samples from 65 patients, 51
pre-treatment and 58 on-treatment, with RECIST response labels (48 PD, 34 SD,
23 PRCR, 4 UNK).

Supplementary files at
`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/`:
`GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz` (and `.raw`, `.rld`
variants). Response labels are in the series matrix at
`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/matrix/GSE91061_series_matrix.txt.gz`.

This cohort is explicitly **not** the spine: 51 pre-treatment samples against
~20,000 genes cannot support honest supervised learning. It is reserved as a
Phase 4 extension for paired pre/on-treatment analysis, where its small size is
acceptable because the paired design supplies statistical power.

### Why this dataset was chosen

Each property below was verified, not assumed, and each maps to a specific
teaching objective:

1. **Long-format clinical data.** The API returns 9,966 patient-attribute rows for
   1,661 patients, requiring a pivot to a usable table.
2. **Everything is a string.** `OS_MONTHS`, `TMB_NONSYNONYMOUS`, and `AGE` all
   arrive as strings; `OS_STATUS` is the composite string `"1:DECEASED"`.
3. **Severe right skew.** TMB has median 5.87 against max 207.5, making log
   transformation necessary rather than decorative.
4. **A genuine batch confound.** Three panels of differing size (341/410/468 genes)
   produce non-comparable raw TMB values. The published correction — ranking TMB
   within cancer type — requires `groupby().transform()`, the highest-value
   intermediate pandas idiom.
5. **Censored survival data.** 829 patients are alive at last contact, not cured.
   Naive handling silently biases results.
6. **Balanced outcome.** 832/829 death/alive supports classification without
   class-imbalance complications masking other lessons.
7. **An honestly weak effect.** TMB is a real but modest predictor. The project
   cannot produce a misleadingly excellent model.

### Known limitations

- Retrospective, non-randomized cohort; no untreated control arm. The analysis
  can establish association only, never causation or treatment benefit.
- Panel heterogeneity is mitigated by within-cancer-type ranking but not
  eliminated.
- No direct immune measurements — no cytokines, cell counts, or expression. The
  immunological reasoning is mechanistic inference from mutation burden. Phase 4
  addresses this.
- Overall survival conflates immunotherapy response with unrelated mortality.

## Architecture

```
immunotherapy-tmb-survival/
├── README.md
├── requirements.txt
├── Makefile
├── .gitignore
├── data/
│   ├── raw/          # API responses; immutable; gitignored
│   ├── interim/
│   └── processed/    # analysis-ready tables
├── src/
│   ├── data/
│   │   ├── fetch.py      # cBioPortal client, caching, retry
│   │   └── clean.py      # pivot, dtype coercion, joins
│   ├── features/
│   │   └── tmb.py        # log transform, within-cancer-type percentile
│   ├── analysis/
│   │   └── survival.py   # Kaplan-Meier, log-rank, Cox
│   ├── models/
│   │   └── train.py      # baseline and regularized models
│   └── viz/
│       └── plots.py      # shared styling
├── notebooks/
│   ├── 01-eda.ipynb
│   ├── 02-survival.ipynb
│   └── 03-modeling.ipynb
├── tests/
├── figures/
└── docs/superpowers/specs/
```

**Central architectural rule:** all logic lives in `src/` and is imported by
notebooks. Notebooks provide narrative and figures; they contain no
transformation logic worth testing. This keeps transformations unit-testable and
lets a reviewer assess code quality without reading notebook cells.

**Module boundaries.** Each module has one responsibility and a documented
interface:

- `fetch.py` — network only. Knows about HTTP and caching; knows nothing about
  what the data means. Writes raw JSON to `data/raw/`, never overwrites, never
  re-downloads when cached.
- `clean.py` — shape only. Consumes raw JSON, produces a tidy patient-level
  DataFrame. Knows nothing about TMB or survival.
- `tmb.py` — the feature transformations, as pure functions on DataFrames.
- `survival.py` — statistics. Consumes the processed table, returns fitted models
  and test statistics.
- `train.py` — modeling. Consumes processed features, returns fitted estimators
  and evaluation metrics.
- `plots.py` — presentation only. Consumes results, produces figures.

Data flows one direction: `fetch → clean → tmb → {survival, train} → plots`. No
module imports a module downstream of it.

## Phases

Each phase is independently shippable. Work abandoned after any phase still
leaves a coherent, presentable project.

### Phase 0 — Acquisition and environment (~1 evening)

Virtual environment, `requirements.txt`, `.gitignore`, repository skeleton, and
`src/data/fetch.py`. Raw API responses cached to `data/raw/`.

**Teaches:** reproducible environments; the raw/processed separation.
**Done when:** a clean clone plus one command produces cached raw data.

### Phase 1 — Wrangling and EDA (~2 weeks)

- Pivot 9,966 long-format rows to a 1,661-row patient table
- Join patient-level (survival, drug, sex) to sample-level (TMB, cancer type,
  panel) data
- Coerce dtypes; parse `OS_STATUS` from `"1:DECEASED"` to boolean
- Plot raw TMB, observe the skew, apply log transform, re-plot
- Compare TMB across the 11 cancer types with `groupby`
- Compute within-cancer-type TMB percentile via `groupby().transform()`
- Benchmark a Python loop against a vectorized numpy operation on the same task

**Teaches:** reshaping, joining, dtypes, split-apply-combine, vectorization,
matplotlib/seaborn.
**Done when:** a notebook contains 6-8 figures with written interpretation.

### Phase 2 — Statistics (~3 weeks)

- Kaplan-Meier curves for high versus low TMB, with log-rank test
- Explicit treatment of censoring: deliberately mishandle it first, observe the
  distortion, then handle it correctly
- Per-cancer-type analysis producing 11 tests, corrected with Benjamini-Hochberg
- Effect sizes with confidence intervals alongside p-values
- Cox proportional hazards adjusting for cancer type and drug class

**Teaches:** survival analysis, hypothesis testing, multiple-testing correction,
scipy and lifelines.
**Done when:** a written statistical report resembling a paper's results section.

### Phase 3 — Machine learning (~3-4 weeks)

**Open design decision, deliberately deferred to the author at implementation
time:** how to define the prediction target. Binary classification ("alive at 12
months") is simpler and maps directly onto coursework but discards patients
censored before 12 months. A survival model uses every patient but is
conceptually harder. This trade-off between statistical honesty and tractability
has no objectively correct resolution and is the author's to make.

- Stratified train/test split, with explanation of why random splitting misleads
- Baseline logistic regression on TMB alone, established before any complex model
- Feature expansion: cancer type, drug class, age, sex, per-gene mutations
  (pivoted from the long-format mutation table)
- Regularization, cross-validation, feature importance
- Honest evaluation against the baseline

**Teaches:** scikit-learn, train/test discipline, overfitting, regularization,
honest evaluation.
**Done when:** a model exists with evaluation against baseline and a written
account of its limitations.

Expected outcome: a modestly performing model. TMB genuinely explains only part
of immunotherapy response. Reporting this accurately is the objective.

### Phase 4 — Extension: Riaz melanoma cohort (optional)

Paired pre/on-treatment RNA-seq from `GSE91061`, addressing a different question:
what changes in the tumor when the drug works? Introduces paired statistical
design and direct immune gene expression.

## Implementation planning

This spec covers roughly ten weeks of work and is **too large for a single
implementation plan**. Each phase gets its own plan, written immediately before
that phase begins rather than all at once.

Planning order:

1. Write an implementation plan for **Phase 0 + Phase 1 together** — they share
   the acquisition-to-clean-data arc and are natural to build in one pass.
2. Execute, review, then plan Phase 2.
3. Execute, review, then plan Phase 3.
4. Phase 4 only if desired after Phase 3 ships.

Planning later phases just-in-time is deliberate: what the author finds difficult
in Phase 1 should change how Phase 2 is scoped, and a plan written ten weeks
early would be based on guesses about their skill level.

## Testing strategy

Approximately 10-15 tests in `tests/`, targeting transformations whose silent
failure would corrupt results:

- The long-to-wide pivot preserves patient count: 1,661 in, 1,661 out
- `OS_STATUS` parsing maps `"1:DECEASED"` to `True`, `"0:LIVING"` to `False`, and
  raises on any unexpected value
- The patient/sample join neither drops nor duplicates patients — an explicit
  regression test, since a silent empty merge is a known failure mode of this
  data shape
- Within-cancer-type percentiles lie in [0, 1] and are approximately uniform
  within each group
- Log transformation handles TMB values of exactly 0.00, which are present

Tests use small synthetic fixtures, not the full dataset, so the suite runs in
seconds without network access.

## Error handling

- **Network.** `fetch.py` retries with backoff, fails loudly on persistent error,
  and never writes a partial cache file. A truncated cache that looks valid is
  worse than no cache.
- **Schema drift.** After loading, assert expected patient count and required
  columns. cBioPortal may change; the pipeline must fail visibly rather than
  produce quietly wrong numbers.
- **Unexpected categorical values.** Parsers raise on unrecognized input rather
  than defaulting. A new `OS_STATUS` value must not silently become `False`.
- **Missing data.** Handled explicitly per column with documented rationale, never
  by a blanket `dropna()`.

## Reproducibility contract

`git clone` → `pip install -r requirements.txt` → `make all` regenerates every
figure and result from scratch.

Makefile targets:

| Target | Effect |
| --- | --- |
| `make data` | Fetch from the API and cache to `data/raw/` (no-op if cached) |
| `make process` | Build analysis-ready tables in `data/processed/` |
| `make figures` | Regenerate all figures into `figures/` |
| `make test` | Run the pytest suite (no network required) |
| `make all` | `data` → `process` → `figures` |

Raw data is gitignored; figures are committed so the README renders on GitHub
without requiring a clone.

## README structure

1. Question and headline finding, with the key figure at the top
2. Data provenance and citation (Samstein et al., *Nature Genetics* 2019, via
   cBioPortal; CC-BY where applicable)
3. Methods summary
4. Results with figures
5. Limitations, stated plainly
6. Reproduction instructions

Section 5 is required and must be specific. Correctly identifying why this
analysis cannot establish causation demonstrates more competence than any
accuracy figure.

## Development practices

- Public repository from the first commit
- Incremental commits throughout the ten weeks; the git history is itself
  evidence of process
- Notebooks committed with outputs cleared where they would otherwise create
  large diffs

## Success criteria

The project succeeds if:

1. A stranger can clone the repository and reproduce every figure
2. The author can explain every transformation in an interview
3. The statistical conclusions are correct, including correct handling of
   censoring and multiple testing
4. Limitations are stated accurately rather than minimized
5. The author is fluent in pandas reshaping, split-apply-combine, and
   vectorization by the end of Phase 1

## Out of scope

- Deep learning
- Raw sequence processing, alignment, or variant calling
- Interactive dashboards or web deployment
- Cohorts beyond MSK-IMPACT and (optionally) Riaz
- Causal inference from this observational data
