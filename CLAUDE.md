# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project context

Retrospective cohort study on peritoneal carcinomatosis patients. The core research question across all analyses is whether surgical PCI (`sPCI`) or pathological PCI (`pPCI`) is the better predictor of overall survival, and whether the two scores can be treated as interchangeable.

Key variables:
- `sPCI` — Peritoneal Cancer Index scored intraoperatively by the surgeon
- `pPCI` — same index scored post-operatively by the pathologist
- `CC` — Completeness of Cytoreduction (binary covariate in Cox models)
- `OS` — Overall Survival (outcome)

## Environment

Managed with [Poetry](https://python-poetry.org/). Always activate the environment before running anything:

```bash
poetry shell        # activate
poetry install      # install/sync dependencies
```

Run a notebook non-interactively:

```bash
poetry run jupyter nbconvert --to notebook --execute notebooks/<name>.ipynb --output notebooks/<name>.ipynb
```

There are no automated tests beyond the `tests/` stub. To verify library functions, run cells interactively in the relevant notebook.

## Architecture

### Analysis pipeline

Notebooks are numbered and should be run in order; each builds on cleaned data from the previous stage:

| Notebook | Purpose |
|----------|---------|
| `01-agreement_all` | Inter-rater agreement (sPCI vs pPCI) — see findings below |
| `02-assoc_difference_*` | OLS, robust (RLM), and quantile regression of sPCI−pPCI on covariates |
| `03-survival_wrangling` | Data cleaning and survival object construction |
| `03-survival_kaplan_meier` / `03-survival_KM_subgroups` | Kaplan-Meier curves |
| `04-cox_standard` / `04-cox_interaction` | Cox PH models with sPCI and pPCI as predictors |
| `05-cox_no_sex` | Sensitivity analysis: Cox model without sex covariate |

Processed datasets used downstream live in `data/`:
- `all_data.csv` (tab-separated, n=454) — full dataset for agreement analysis
- `model_A_data.csv` / `model_B_data.csv` — subsets used in Cox models
- `GPT_processed_survival_data*.csv` — GPT-assisted data cleaning output

### Custom library (`src/survival/agreement.py`)

All non-standard statistical functions are centralised here and imported into notebooks. The library is installed as an editable package via Poetry (`packages = [{include = "survival", from = "src"}]`).

Key functions:

| Function | Description |
|----------|-------------|
| `paired_permutation_test` | Sign-flip permutation test on mean difference |
| `paired_permutation_test_median` | Same, for median difference |
| `paired_bootstrap_ci` / `paired_bootstrap_ci_bca` | Percentile and BCa bootstrap CIs |
| `sign_test` | Exact binomial sign test |
| `bland_altman_analysis` | Classical Bland-Altman plot + LoA |
| `bland_altman_analysis_robust` | Quantile-based LoA + proportional bias regression |
| `bland_altman_regression_loa` | Regression-based LoA for proportional + heteroscedastic bias |
| `intraclass_correlation` | ICC(A,1) via `pingouin` |
| `icc_bias_corrected` | ICC after removing median fixed offset |

All functions follow the convention `(method1, method2, ...)` where `method1=sPCI`, `method2=pPCI`.

Plots are saved to `plots/` at 600 dpi. Reports are Markdown/PDF in `reports/`.

## Key findings — Agreement analysis (notebook 01)

sPCI and pPCI **cannot be treated as interchangeable**.

### Fixed bias
sPCI exceeds pPCI by ~5–7 units on average (mean diff=6.94, median diff=5.0). Sign test: 399/429 pairs have sPCI > pPCI (p≈2×10⁻⁸³).

### Proportional bias
Regression of differences on means: `d = 0.913 + 0.486·mean` (r=0.514, p=6.1×10⁻³²). The gap widens at higher disease burden.

### Regression-based limits of agreement (Bland & Altman 1999)
Residuals are also heteroscedastic (`|residuals| ~ mean`: slope=0.177, p=1.1×10⁻¹¹), so the LoA are not constant. Fixed LoA misrepresent agreement at both extremes.

| Mean PCI | Bias  | LoA lower | LoA upper |
|----------|-------|-----------|-----------|
| 5        |  3.34 |    −3.81  |    10.49  |
| 10       |  5.78 |    −3.54  |    15.09  |
| 15       |  8.21 |    −3.28  |    19.69  |
| 20       | 10.64 |    −3.01  |    24.29  |
| 25       | 13.07 |    −2.75  |    28.89  |

The upper LoA grows from ~10.5 to ~28.9; the lower LoA stays near zero because pPCI almost never exceeds sPCI.

### ICC
- Raw ICC(A,1) = 0.46 — poor absolute agreement
- After removing median bias (δ=5): ICC = 0.62 [0.54, 0.68] — moderate consistency in relative rank, but proportional bias is not corrected

## Cox model findings (notebooks 04–05)

n=263 patients, 100 death events, median follow-up 13 months (truncated at 57 months). Both sPCI and pPCI are independent predictors of OS, but models 04 vs 05 compare their relative predictive value. The comparison report is in `reports/05-vs-04_cox_comparison.md/pdf`.
