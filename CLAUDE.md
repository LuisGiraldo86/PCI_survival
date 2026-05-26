# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project context

Retrospective cohort study on peritoneal carcinomatosis patients. The core research question across all analyses is whether surgical PCI (`sPCI`) or pathological PCI (`pPCI`) is the better predictor of overall survival, and whether the two scores can be treated as interchangeable.

Key variables:
- `sPCI` — Peritoneal Cancer Index scored intraoperatively by the surgeon
- `pPCI` — same index scored post-operatively by the pathologist
- `CC` — Completeness of Cytoreduction (0 = complete, 1 = incomplete)
- `OS` — Overall Survival (outcome)
- `Datum Rezidiv` — despite the name ("recurrence date"), this column records the **last date the patient was seen alive** (last clinical contact). It is used as the censoring date for patients who did not die. Censoring is non-informative.

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
| `02.1-assoc_difference_ols` | OLS regression of sPCI−pPCI on covariates (baseline) |
| `02.2-assoc_difference_ols_robust` | Robust regression (RLM, IRLS Huber-T) — preferred inference |
| `02.3-assoc_difference_quantile` | Quantile regression across q = 0.25, 0.50, 0.75, 0.90 |
| `03.1-survival_wrangling` | Data cleaning and survival object construction |
| `03.2-survival-kaplan_meier` | Overall KM curve, truncation rule, RMST |
| `03.3-survival_KM_subgroups` | KM curves by subgroup (sex, CC, tumour type, age) |
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

## Key findings — discordance analysis (notebooks 02.1, 02.2, 02.3)

**N = 409** after excluding tumour type 8 and CC > 1. Outcome: `raw_diff` = sPCI − pPCI. Reference: tumour type 1. Both notebooks address the same model; OLS is the baseline, RLM is the preferred inference.

### OLS (baseline)
R² = 0.373, Adj R² = 0.356. Assumptions violated: Shapiro-Wilk W = 0.980 (p < 0.001), Breusch-Pagan p < 0.001. OLS standard errors are unreliable; p-values may be anti-conservative.

### RLM — preferred model (IRLS, Huber-T, bootstrap CIs)
Weighted R² = 0.383. Results consistent across HuberT, Bisquare, and AndrewWave norms.

| Predictor | OLS coef | RLM coef | Bootstrap 95 % CI | p (RLM) |
|---|---|---|---|---|
| Tumour type 5 | +6.79 | +5.88 | [+4.10, +7.98] | < 0.001 |
| CC (incomplete) | +6.02 | +5.53 | [+4.08, +7.15] | < 0.001 |
| Tumour type 7 | +4.45 | +4.18 | [+1.74, +6.83] | < 0.001 |
| Thermoablation | −1.39 | −1.62 | [−3.13, −0.08] | 0.032 |

Sex, Age, PreOP CTx, and tumour types 2, 3, 4, 6 are not significant in either model.

**Thermoablation is the weakest finding** — OLS p = 0.047 under violated homoscedasticity, RLM bootstrap CI barely excludes zero. Interpret cautiously; likely reflects patient selection rather than a causal effect.

### Open methodological questions
- CC is downstream of sPCI (collider/mediator risk); a model without CC should be run for comparison.
- Tumour type effects may partly reflect proportional bias (`d ≈ 0.913 + 0.486·mean`); adjusting for sPCI or mean PCI would isolate independent type effects.
- Consider collapsing 02.1 and 02.2 into one: OLS as baseline, RLM as the main model, with no loss of content.

## Key findings — quantile regression (notebook 02.3-assoc_difference_quantile)

Same N=409, same predictors. Quantile regression at q = 0.25, 0.50, 0.75, 0.90. Two models: full (with Age) and reduced (Age dropped — pseudo-R² loss of 0.0003, negligible). Reduced model is preferred.

### Median model (q = 0.5) — reduced model, bootstrap CIs
Pseudo-R² = 0.189.

| Predictor | Coef | Bootstrap 95 % CI | Significant |
|---|---|---|---|
| CC (incomplete) | +5.00 | [+3.00, +7.00] | Yes |
| Tumour type 5 | +6.00 | [+3.00, +8.00] | Yes |
| Tumour type 7 | +5.00 | [+1.00, +8.00] | Yes |
| Thermoablation | −2.00 | [−4.00, −0.48] | Yes |
| PreOP CTx | −2.00 | [−3.00, +0.00] | Ambiguous — do not report as confirmed |

### Key new insight: CC and tumour type 5 effects are asymmetric
| Predictor | q = 0.25 | q = 0.50 | q = 0.75 | q = 0.90 | Asymmetric? |
|---|---|---|---|---|---|
| CC | +3.6*** | +5.2*** | +9.0*** | +9.1*** | Yes — CI [+1.9, +7.9] |
| Tumour type 5 | sig | sig | sig | sig | Yes — CI [+0.8, +7.9] |
| Thermoablation | −1.4* | −2.1** | NS | NS | No |

CC effect nearly doubles from lower to upper quantile — incomplete cytoreduction is disproportionately associated with discordance in patients already at the high end. OLS/RLM's single coefficient masks this gradient. Thermoablation is significant only at the lower two quantiles; the protective association disappears in high-discordance patients.

### Technical limitation
`raw_diff` is integer-valued, causing ties that break the quantile regression sparsity estimator and collapse bootstrap CIs to round numbers (e.g., PreOP CTx CI [−3.00, +0.00]). Repeated `IterationLimitWarning` during bootstrap is a direct consequence.

### Planned next steps
1. ~~**Jitter outcome** — add `U(−0.5, 0.5)` noise to break ties; refit and check stability across seeds (Machado & Silva 2005)~~ → done in 02.4
2. ~~**Model without CC** — quantify collider/mediator bias on tumour type coefficients~~ → done in 02.4
3. ~~**Add mean PCI covariate** — isolate type effects from proportional bias mechanism~~ → done in 02.4
4. ~~**Drop PreOP CTx** from final model — not confirmed across any of the three modelling approaches~~ → revised in 02.4 (see below)
5. **Test CC × thermoablation interaction** — thermoablation selection may be confounded with CC status
6. **Consider negative binomial on abs\_diff** if jittering proves unstable

## Key findings — sensitivity analyses (notebook 02.4)

### Action 1 — Jitter `raw_diff`

Jittering with U(−0.5, +0.5) is stable across 5 seeds (max range per predictor: type 3 = 0.78, CC = 0.56, type 5 = 0.26). Bootstrap CIs are no longer degenerate.

Jittered bootstrap CIs (q=0.5, seed=42, n=2000):

| Predictor | Coef | 95 % CI |
|---|---|---|
| CC | +5.06 | [+3.46, +6.86] |
| Tumour type 5 | +5.79 | [+3.36, +7.88] |
| Tumour type 7 | +5.23 | [+1.04, +7.67] |
| Thermoablation | −2.05 | [−3.90, −0.51] |
| PreOP_CTx | −1.63 | [−2.96, +0.04] |

### Action 4 — PreOP_CTx: revised verdict

The degenerate CI in 02.3 was a numerical artefact. With jittering, PreOP_CTx is consistently p=0.016–0.040 at q=0.5 across all 5 seeds. Bootstrap CI [−2.96, +0.04] barely includes zero. NS at q=0.25, q=0.75, q=0.90. **Borderline median-only predictor** — not reportable as confirmed, but prior dismissal was premature. A sensitivity model including it is warranted.

### Action 2 — Model without CC

R² drops from 0.371 → 0.217 (OLS) when CC is removed — CC carries large independent variance.

| Predictor | With CC | No CC | Δ OLS | Δ RLM | Δ QR |
|---|---|---|---|---|---|
| Tumour type 5 | +6.67 | +7.01 | +0.34 | −0.06 | −0.15 |
| Tumour type 7 | +4.48 | +5.43 | +0.95 | +1.02 | +1.93 |
| Thermoablation | −1.42 | −2.53 | −1.11 | −1.11 | +0.09 |

Tumour type 5 is essentially unconfounded by CC. Tumour type 7 inflates ~1–2 units without CC — mild confounding, effect survives. Thermoablation strengthens notably without CC, consistent with selection bias (thermoablation patients disproportionately achieve CC0).

### Action 3 — Add mean_pci

`mean_pci = (sPCI + pPCI) / 2` is strongly significant across all models (OLS: +0.314***, RLM: +0.347***, QR: +0.315***) and grows toward the upper tail in QR (+0.215 at q=0.25 → +0.475 at q=0.75). R²/pseudo-R² gain of +0.07 in both OLS and QR.

| Predictor | Without mean_pci | With mean_pci | Δ (QR) |
|---|---|---|---|
| CC | +4.94 | +2.64 | −2.30 |
| Tumour type 5 | +5.41 | +4.90 | −0.50 |
| Tumour type 7 | +4.76 | +2.96 | −1.81 |
| Thermoablation | −2.19 | −1.35 | +0.84 |

**Tumour type 5 is the most robust finding** — survives CC removal and mean_pci adjustment with < 10% attenuation. CC halves after mean_pci adjustment but remains significant. Tumour type 7 attenuates ~38% in QR — partly a burden proxy, but survives.

### Overall picture

sPCI−pPCI discordance is driven by three independent mechanisms:
1. **Disease burden** (mean_pci, +0.31 per unit) — proportional bias confirmed in regression
2. **Tumour biology** (type 5 most robust; type 7 with mild CC confounding)
3. **Surgical outcome** (CC, asymmetric — effect nearly doubles from lower to upper quantile)

Thermoablation remains fragile (median-only, partly CC-confounded). The CC asymmetry (q=0.25: +3.6 → q=0.90: +9.1) is the most clinically meaningful result: discordance in CC1 patients is disproportionately worse in those already at the high end.

## Cox model findings (notebooks 04–05)

n=263 patients, 100 death events (38%). Median OS = 35 months. Median follow-up = 13 months (median observation time across all patients including early-censored). Truncated at τ = 57 months (first point where < 10% at risk). RMST = 34.59 months (bootstrap 95% CI: [31.63, 37.55]). Both sPCI and pPCI are independent predictors of OS, but models 04 vs 05 compare their relative predictive value. The comparison report is in `reports/05-vs-04_cox_comparison.md/pdf`.
