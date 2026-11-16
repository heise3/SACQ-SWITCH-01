# SACQ-SWITCH-01

Analysis code, frozen results and figure source data accompanying **Humoral and granulocyte programs diverge across lupus activity states, treatment and renal compartments** by Peng Liang, Jie Wu, BaoLong Li, XiangYan Wu and Xiaona Zhu.

This release combines component-level lupus activity scores, within-patient state contrasts, prospective models, adult treatment trajectories, donor-level single-cell profiles and renal-compartment comparisons. The inferential unit is the patient or donor. Cell counts are descriptive. Results and model estimates are unchanged from the completed scientific analysis; version 1.3.2 revises figure layout and citation/submission preparation only.

## Contents

- `code/task1.py` to `task4.py`: original analysis implementations.
- `code/common.py`, `phase0_parsers.py`: portable public-input adapter and fixed module definitions. The parser declarations are separated from the historical runner; no private workspace is required.
- `results/`: frozen results and derived display tables, including public study sample/donor labels. No direct identifiers are present.
- `figures/`: four main and two supplementary figures, each editable SVG, vector PDF and 300 dpi PNG.
- `tables/`: machine-readable main tables.
- `references/`: verified article metadata and CSL records.
- `documentation/EXTERNAL_INPUTS.json`: official input URLs, sizes and SHA-256 values.
- `documentation/FIGURE_SOURCE_DATA_MAP_v1.3.json`: panel-level source and inference mapping. Historical `supplement/` paths in that map refer to `documentation/` in this release; the geometry is `KM_GEOMETRY.json`.

## Reproduce figures without refitting

Create a Python environment from `documentation/requirements.lock.txt`, then run from the checkout root:

```sh
python code/verify_release.py
python code/render_figures.py
```

New drawings are written to `rendered_figures/`; shipped figures and frozen results are retained. Kaplan–Meier geometry is reused from the scientific analysis. Available-unit means and descriptive pointwise t intervals use the original display formula. No analysis model or inferential resampling is executed by the renderer.

## Full reanalysis

Full reanalysis additionally needs the original expression inputs (not redistributed). `python code/fetch_inputs.py` retrieves only the documented public GEO and CELLxGENE assets and checks analyzed-file hashes. Download into an isolated checkout on storage with sufficient capacity. Current repository contents are not an expression-data bundle. Upstream file changes cause a hash failure and must be investigated before any reanalysis.

For an intentional fresh analysis, use a separate copy, preserve this release, create `work/`, and run `python code/task1.py`, then tasks 2, 3 and 4 in order. Task 1 reuses `results/task1_checkpoint.json` when present; preserve that checkpoint and remove it only in the separate reanalysis copy if new bootstrap computation is intended. Native numerical-library and software-version differences may affect newly fitted results. Full raw-input replay was not run during this editorial release. The portable adapter was statically checked; supplied-result figure replay was checked independently.

This exploratory study does not identify causal mediation, establish selective absolute neutrophil depletion, or validate a clinical risk predictor. The adult induction cohort is not an independent adult SACQ incident-activity validation cohort. Nominal and multiplicity-adjusted findings are distinguished in the results.

## Citation and rights

Cite the versioned Zenodo archival record and the original cohort publications. The archived v1.3.2 release is published at https://doi.org/10.5281/zenodo.24294834 (Git commit `18d984cf0c19e763298c666233e0411496b3ad4c`). The GitHub source repository is https://github.com/heise3/SACQ-SWITCH-01.

Copyright 2026 the authors. Source data, repository annotations, fonts and third-party software retain their own terms. Author-created code and derived outputs are licensed under CC BY 4.0, matching the published Zenodo record; no license is conferred on third-party inputs by this repository. The manuscript and submission correspondence are not included.
