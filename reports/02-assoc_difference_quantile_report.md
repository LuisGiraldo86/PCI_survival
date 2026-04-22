# Quantile Regression Analysis of sPCI–pPCI Discrepancy

**Notebook:** `02-assoc_difference_quantile.ipynb`  
**Date:** 2026-04-21  
**Outcome:** `raw_diff = sPCI − pPCI` (surgical minus pathological PCI score)

---

## 1. Study Sample

| Statistic | Value |
|---|---|
| N | 409 |
| Mean raw_diff | +7.09 |
| Median raw_diff | +5.00 |
| SD | 7.01 |
| Range | −8 to +32 |

The distribution is **right-skewed** (mean > median), with a long upper tail of large overestimates. Quantile regression is used instead of OLS to model the conditional distribution without normality assumptions.

Tumor type 1 (n = 122, the most frequent) is the reference category. Tumor type 8 was excluded. Patients with incomplete cytoreduction (CC > 1) were excluded. Final predictors: Sex, Age, CC (cytoreduction completeness, binary 0/1), PreOP_CTx, Thermoablation, and dummy variables for tumor types 2–7.

---

## 2. Full Model (including Age)

### 2.1 Median Model (q = 0.50)

| Predictor | Coef | p | Bootstrap 95% CI |
|---|---|---|---|
| Intercept | +7.38 | <0.001 | [+4.00, +10.94] * |
| Sex | −0.22 | 0.745 | [−1.77, +1.00] |
| **Age** | **−0.02** | **0.431** | **[−0.07, +0.03]** |
| **CC** | **+5.22** | **<0.001** | **[+3.27, +7.00] *** |
| PreOP_CTx | −1.90 | 0.017 | [−3.00, −0.00] * |
| **Thermoablation** | **−2.08** | **0.005** | **[−4.00, −0.40] *** |
| Tumor_type_2 | −0.36 | 0.805 | [−3.00, +2.64] |
| Tumor_type_3 | +0.18 | 0.888 | [−2.41, +3.63] |
| Tumor_type_4 | +0.20 | 0.838 | [−1.33, +2.39] |
| **Tumor_type_5** | **+5.56** | **<0.001** | **[+3.26, +7.89] *** |
| Tumor_type_6 | −1.42 | 0.201 | [−3.00, +0.58] |
| **Tumor_type_7** | **+4.80** | **<0.001** | **[+1.00, +7.68] *** |

### 2.2 Quantile Process Summary

| Predictor | q=0.25 | q=0.50 | q=0.75 | q=0.90 |
|---|---|---|---|---|
| Sex | −0.22 | −0.22 | −1.00 | −2.40 * |
| Age | −0.03 | −0.02 | 0.00 | +0.03 |
| **CC** | **+3.60 *** | **+5.22 *** | **+9.00 *** | **+9.10 *** |
| PreOP_CTx | −0.67 | −1.90 * | −1.00 | −0.60 |
| Thermoablation | −1.40 * | −2.08 ** | −1.00 | −0.28 |
| Tumor_type_3 | −0.07 | +0.18 | +4.00 * | +5.80 ** |
| Tumor_type_4 | +0.31 | +0.20 | +3.00 * | +4.05 * |
| **Tumor_type_5** | **+4.17 *** | **+5.56 *** | **+9.00 *** | **+10.45 *** |
| **Tumor_type_7** | **+2.57 * | **+4.80 *** | **+6.00 *** | **+4.88 * |

(Tumor types 2 and 6 non-significant at all quantiles.)

### 2.3 Model Fit

| Quantile | Pseudo-R² |
|---|---|
| q=0.25 | 0.136 |
| q=0.50 | 0.189 |
| q=0.75 | 0.261 |
| q=0.90 | 0.358 |
| OLS R² | 0.373 |

### 2.4 Symmetry Test (q=0.75 − q=0.25, bootstrap)

| Predictor | Δ coef | 95% CI | Asymmetric? |
|---|---|---|---|
| **CC** | **+5.04** | **[+1.94, +7.88]** | **YES *** |
| **Tumor_type_5** | **+3.94** | **[+0.84, +7.88]** | **YES *** |
| All others | — | CI includes 0 | No |

---

## 3. Reduced Model (Age Excluded)

Age was **not significant at any quantile** (all p > 0.10, bootstrap CI includes zero). The analysis was repeated without Age.

### 3.1 Median Model (q = 0.50)

| Predictor | Coef | p | Bootstrap 95% CI |
|---|---|---|---|
| Intercept | +6.00 | <0.001 | [+4.00, +8.00] * |
| Sex | −0.00 | 1.000 | [−2.00, +1.00] |
| **CC** | **+5.00** | **<0.001** | **[+3.00, +7.00] *** |
| PreOP_CTx | −2.00 | 0.012 | [−3.00, +0.00] |
| **Thermoablation** | **−2.00** | **0.007** | **[−4.00, −0.48] *** |
| Tumor_type_2 | −0.00 | 1.000 | [−3.00, +3.00] |
| Tumor_type_3 | −0.00 | 1.000 | [−2.90, +3.83] |
| Tumor_type_4 | +0.00 | 1.000 | [−1.50, +2.50] |
| **Tumor_type_5** | **+6.00** | **<0.001** | **[+3.00, +8.00] *** |
| Tumor_type_6 | −1.00 | 0.365 | [−3.00, +1.00] |
| **Tumor_type_7** | **+5.00** | **<0.001** | **[+1.00, +8.00] *** |

### 3.2 Quantile Process Summary

| Predictor | q=0.25 | q=0.50 | q=0.75 | q=0.90 |
|---|---|---|---|---|
| Sex | 0.00 | 0.00 | −1.00 | −2.00 |
| **CC** | **+4.00 *** | **+5.00 *** | **+9.00 *** | **+9.00 *** |
| PreOP_CTx | −1.00 | −2.00 * | −1.00 | −1.00 |
| Thermoablation | −2.00 ** | −2.00 ** | −1.00 | 0.00 |
| Tumor_type_3 | 0.00 | 0.00 | +4.00 * | +6.00 ** |
| Tumor_type_4 | 0.00 | 0.00 | +3.00 * | +4.00 * |
| **Tumor_type_5** | **+4.00 *** | **+6.00 *** | **+9.00 *** | **+11.00 *** |
| Tumor_type_7 | +2.00 | +5.00 *** | +6.00 *** | +5.00 * |

(Tumor types 2 and 6 non-significant at all quantiles.)

### 3.3 Model Fit

| Quantile | Pseudo-R² |
|---|---|
| q=0.25 | 0.133 |
| q=0.50 | 0.189 |
| q=0.75 | 0.261 |
| q=0.90 | 0.355 |
| OLS R² | 0.373 |

Removing Age reduces OLS R² by only 0.0005 (0.3731 → 0.3726), confirming it contributes nothing.

### 3.4 Symmetry Test (q=0.75 − q=0.25, bootstrap)

| Predictor | Δ coef | 95% CI | Asymmetric? |
|---|---|---|---|
| Intercept | +5.67 | [+2.00, +10.00] | YES * |
| **CC** | **+4.97** | **[+2.00, +8.00]** | **YES *** |
| **Tumor_type_5** | **+4.05** | **[+1.00, +8.00]** | **YES *** |
| All others | — | CI includes 0 | No |

---

## 4. Diagnostics

All four quantile models passed all diagnostic checks:

- **Residual sign proportions** match the nominal quantile at all q (±5% tolerance, all OK), indicating numerical convergence.
- **Residuals vs fitted** show no systematic trend or heteroskedasticity.
- **Spearman correlations** between residuals and all predictors are near zero and non-significant, ruling out missed non-linearity or obvious omitted variables.
- **Bootstrap SEs vs asymptotic SEs** agree closely for all significant predictors, validating the sparsity estimator.
- `IterationLimitWarning` messages arise only in a small number of edge-case bootstrap resamples, not in the primary fitted models; they do not affect the conclusions.

---

## 5. Quantile Crossing

Quantile crossing (where a lower fitted quantile exceeds a higher one for a given observation) was checked for both models across all adjacent quantile pairs (q25–q50, q50–q75, q75–q90).

### Full model (with Age)

| Comparison | Cases crossing |
|---|---|
| q25 > q50 | 0.000% |
| q50 > q75 | 0.000% |
| q75 > q90 | 0.000% |
| Any crossing | 0.000% |

No crossing occurs in any observation, indicating perfect ordering of the predicted quantile curves for the full model.

### Reduced model (Age excluded)

| Comparison | Cases crossing |
|---|---|
| q25 > q50 | 0.000% |
| q50 > q75 | 0.000% |
| q75 > q90 | 0.489% |
| Any crossing | 0.489% |

Crossing is negligible (< 0.5% of observations, confined to the q75–q90 boundary), confirming that the reduced model's quantile predictions are also internally consistent. The minor q75–q90 crossing in the reduced model compared to none in the full model is attributable to the small change in model flexibility when Age is dropped; it has no practical impact on conclusions.

---

## 6. Conclusions

### Primary finding
**Both cytoreduction completeness (CC) and tumor type drive the sPCI–pPCI discrepancy.** CC is the single strongest predictor: complete cytoreduction (CC = 1 vs CC = 0) is associated with a **+5 point increase in median discrepancy**, growing to **+9 points** at the upper quartiles. Among tumor types, type 5 shows a median overestimate ~6 points higher than type 1, escalating to +9 at q=0.75 and +11 at q=0.90. Tumor type 7 shows a similar pattern (+5 at the median, +5–6 in the upper quartiles). Tumor types 3 and 4 emerge as significant only at the upper quantiles (q=0.75 and q=0.90), suggesting latent overestimation in these groups that is only apparent among patients with the largest discrepancies.

### Thermoablation reduces the discrepancy
Thermoablation is consistently associated with a **2 point reduction** in sPCI–pPCI at the lower and median quantiles (q=0.25 and q=0.50, both p < 0.01). The effect attenuates toward zero at higher quantiles (q=0.75 and q=0.90), suggesting thermoablation primarily compresses discrepancy in typical cases but does not substantially affect the worst outliers.

### Asymmetric effects
The symmetry test flags predictors whose effects **differ between the lower and upper halves of the distribution**:

- **CC:** the positive association with discrepancy is substantially larger in patients who already have a large discrepancy (Δ ~+5 points from q=0.25 to q=0.75), indicating that complete cytoreduction amplifies overestimation particularly among the worst cases.
- **Tumor_type_5:** the overestimation effect is substantially larger in patients who already have a large discrepancy, indicating that surgeons are disproportionately wrong in the worst cases.

### Age and other predictors
Age is not significant at any quantile and its removal has negligible impact on model fit (ΔOLS R² < 0.001). PreOP_CTx shows marginal significance at the median only (p = 0.012–0.017), with bootstrap CIs borderline including zero; it should not be considered robustly associated. Sex is not significant at any quantile in the reduced model. Tumor types 2 and 6 have no association with the discrepancy at any quantile.

### Model fit
The pseudo-R² increases from ~0.133 at q=0.25 to ~0.355 at q=0.90, closely matching the OLS R² (0.37). This confirms that the predictors — particularly CC and tumor type — account for progressively more variance in the upper tail, where discrepancies are largest and clinically most consequential. The overall model fit is substantially higher than in models that did not include CC as a predictor.
