# Predictors of sPCI–pPCI Discordance: OLS and Robust Regression Analysis

**Date:** April 21, 2026  
**Sources:** `notebooks/02-assoc_difference_ols.ipynb`, `notebooks/02-assoc_difference_ols_robust.ipynb`

---

## 1. Background and Objective

A systematic positive bias between intraoperative surgeon PCI (`sPCI`) and preoperative/pathological PCI (`pPCI`) has been established (mean sPCI − pPCI ≈ 6.9, see report `01-agreement_all_report.md`). This analysis investigates **which patient and clinical characteristics explain the magnitude of this discordance**, using ordinary least squares (OLS) regression.

---

## 2. Study Sample

After excluding tumour type 8 and restricting to CC 0/1, the analytic sample comprised **N = 409** patients.  
Tumour type 1 (n = 122) was set as the reference category.

**Predictors included:**

| Variable | Type | Notes |
|---|---|---|
| Sex | Binary | 0 = female, 1 = male |
| Age | Continuous | Years |
| CC | Binary | 0 = complete cytoreduction, 1 = incomplete |
| PreOP CTx | Binary | Any preoperative chemotherapy (≥1 cycle) |
| Thermoablation | Binary | Intraoperative thermoablation performed |
| Tumour type | Nominal | Types 2–7 vs. reference type 1 |

---

## 3. Outcome Variables

Two outcomes were constructed:

- **`raw_diff`** = sPCI − pPCI (signed discordance)
- **`abs_diff`** = |sPCI − pPCI| (magnitude of discordance)

| Statistic | raw\_diff | abs\_diff |
|---|---|---|
| Mean | 7.09 | 7.39 |
| SD | 7.01 | 6.69 |
| Median | 5.0 | 5.0 |
| Min | −8 | 0 |
| Max | 32 | 32 |

The positive mean and median of `raw_diff` confirm that surgeons consistently find **more** peritoneal disease intraoperatively than predicted preoperatively (median discordance of 5 PCI points). Both distributions depart significantly from normality (Shapiro-Wilk: p < 0.001 for both).

---

## 4. OLS Regression Results (standard)

**Outcome: `raw_diff` (sPCI − pPCI)**

| | F-statistic | df | p-value |
|---|---|---|---|
| Overall model | 21.48 | (11, 397) | 3.0 × 10⁻³⁴ |

**Variance explained:** R² = 0.373, Adjusted R² = 0.356

The model is highly significant overall and explains approximately **37% of the variance** in discordance.

### Regression Coefficients

| Predictor | Coefficient | SE | t | p | 95% CI |
|---|---|---|---|---|---|
| Intercept | +6.23 | 1.69 | 3.69 | 0.0003 | [2.90, 9.55] |
| Sex (male) | −1.15 | 0.64 | −1.79 | 0.074 | [−2.40, 0.11] |
| Age | −0.01 | 0.02 | −0.57 | 0.568 | [−0.06, 0.03] |
| **CC (incomplete)** | **+6.02** | 0.61 | 9.83 | **< 0.001** | [4.81, 7.22] |
| PreOP CTx | −0.85 | 0.75 | −1.13 | 0.260 | [−2.32, 0.63] |
| **Thermoablation** | **−1.39** | 0.70 | −1.99 | **0.047** | [−2.76, −0.02] |
| Tumour type 2 | +0.05 | 1.36 | 0.04 | 0.972 | [−2.63, 2.73] |
| Tumour type 3 | +1.17 | 1.21 | 0.97 | 0.335 | [−1.21, 3.54] |
| Tumour type 4 | +1.64 | 0.92 | 1.78 | 0.076 | [−0.17, 3.45] |
| **Tumour type 5** | **+6.79** | 0.82 | 8.23 | **< 0.001** | [5.17, 8.41] |
| Tumour type 6 | −0.48 | 1.05 | −0.46 | 0.647 | [−2.54, 1.58] |
| **Tumour type 7** | **+4.45** | 1.23 | 3.61 | **0.0003** | [2.03, 6.87] |

### Key Findings

- **CC (incomplete cytoreduction)** is the strongest predictor: patients with incomplete cytoreduction (CC1) show on average **+6.0 points** more discordance than those with complete resection (CC0, p < 0.001).
- **Tumour type 5** shows **+6.8 points** more discordance than type 1 (p < 0.001).
- **Tumour type 7** shows **+4.5 points** more discordance than type 1 (p < 0.001).
- **Thermoablation** is associated with **1.4 points less** discordance (p = 0.047), possibly reflecting more localised and predictable disease in patients selected for this procedure.
- Sex, Age, PreOP CTx, and tumour types 2, 3, 4, and 6 are **not significant** predictors.

---

## 5. Robust Regression Results (Huber M-estimator)

Given the OLS violations (non-normality, heteroskedasticity), the model was re-estimated using an **Iteratively Reweighted Least Squares (IRLS) Huber M-estimator** (default HuberT, tuning constant 1.345), which downweights observations with large residuals. Bootstrap confidence intervals (2 000 resamples, seed 42) are reported alongside asymptotic standard errors.

**Outcome: `raw_diff` (sPCI − pPCI)**

**Weighted R²:** 0.383  |  **Scale (MAD):** 5.06  |  **Durbin-Watson:** 1.99

### RLM Coefficients

| Predictor | Coefficient | Asym. SE | p | 95% CI (bootstrap) |
|---|---|---|---|---|
| Intercept | +7.08 | 1.65 | < 0.001 | [+3.77, +10.55] |
| Sex (male) | −0.82 | 0.63 | 0.191 | [−2.06, +0.40] |
| Age | −0.02 | 0.02 | 0.329 | [−0.07, +0.02] |
| **CC (incomplete)** | **+5.53** | 0.60 | **< 0.001** | [+4.08, +7.15] |
| PreOP CTx | −1.15 | 0.74 | 0.119 | [−2.58, +0.22] |
| **Thermoablation** | **−1.62** | 0.68 | **0.018** | [−3.13, −0.08] |
| Tumour type 2 | −0.05 | 1.33 | 0.970 | [−2.28, +2.35] |
| Tumour type 3 | +0.56 | 1.18 | 0.639 | [−1.84, +3.18] |
| Tumour type 4 | +1.25 | 0.90 | 0.166 | [−0.33, +2.93] |
| **Tumour type 5** | **+5.88** | 0.81 | **< 0.001** | [+4.10, +7.98] |
| Tumour type 6 | −0.92 | 1.03 | 0.370 | [−2.57, +0.73] |
| **Tumour type 7** | **+4.18** | 1.21 | **< 0.001** | [+1.74, +6.83] |

### Key Findings (RLM)

- Results are consistent with OLS in direction and significance for all predictors.
- **CC (incomplete)** remains the strongest predictor: +5.5 points discordance vs. OLS +6.0; bootstrap CI excludes zero.
- **Tumour types 5 and 7** retain strong significance (+5.9 and +4.2 points respectively).
- **Thermoablation** −1.6 points (p = 0.018); bootstrap CI just excludes zero [−3.13, −0.08].
- Sex, Age, PreOP CTx no longer approach significance under robust estimation.
- 20 high-leverage observations (h > 2p/n = 0.059); Cook's D never exceeds 1 (max = 0.042), indicating no individually dominant observations.

---

## 6. Residual Diagnostics (OLS)

| Test | Statistic | p-value | Conclusion |
|---|---|---|---|
| Shapiro-Wilk (residuals) | W = 0.980 | < 0.001 | Residuals non-normal |
| Jarque-Bera | 20.16 | < 0.001 | Confirms non-normality |
| Breusch-Pagan (LM) | 64.88 | < 0.001 | Heteroskedasticity present |
| Breusch-Pagan (F) | 6.80 | < 0.001 | Heteroskedasticity present |
| Skewness | 0.48 | — | Mild right skew |
| Kurtosis | 3.54 | — | Slightly leptokurtic |

Both OLS assumptions of **normality** and **homoskedasticity** of residuals are violated. The reported standard errors and p-values should be interpreted with caution; they are likely anti-conservative (i.e., p-values may be overstated).

---

## 7. Limitations and Recommendations

1. ~~**Heteroskedasticity-robust standard errors** (HC3)~~ **Huber M-estimator with bootstrap CIs** has been applied (see Section 5) and yields the preferred inference.
2. **Quantile regression** on the median would be more robust given the skewed, heteroskedastic outcome and would not rely on distributional assumptions.
3. The signed outcome (`raw_diff`) assumes direction matters. If only magnitude is of interest, rerunning the model on `abs_diff` should be considered.
4. Tumour types with small cell counts (types 2, 6, 7) warrant cautious interpretation due to limited statistical power.
5. The negative coefficient for **Thermoablation** likely reflects **selection bias** rather than a causal effect — patients undergoing thermoablation may have more resectable, well-demarcated disease that is easier to assess preoperatively.

---

## 8. Summary

Both OLS (R² = 0.37) and Huber robust regression (weighted R² = 0.38) yield consistent conclusions across all predictors. **CC status** is the dominant driver (~+5.5 to +6.0 points for incomplete vs. complete cytoreduction). **Tumour types 5 and 7** carry substantially more intraoperative surprise than type 1. **Thermoablation** is associated with slightly less discordance (−1.4 to −1.6 points). Sex, age, and preoperative chemotherapy are not independently significant. The robust model confirms these findings are not driven by outliers (max Cook's D = 0.042, Durbin-Watson ≈ 2.0). OLS assumption violations (non-normality, heteroskedasticity) do not change the substantive interpretation, but bootstrap confidence intervals from the robust model are preferred for inference.
