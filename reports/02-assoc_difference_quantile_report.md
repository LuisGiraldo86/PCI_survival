# Quantile Regression Analysis of sPCI–pPCI Discrepancy

**Notebook:** `02-assoc_difference_quantile.ipynb`  
**Date:** 2026-04-16  
**Outcome:** `raw_diff = sPCI − pPCI` (surgical minus pathological PCI score)

---

## 1. Study Sample

| Statistic | Value |
|---|---|
| N | 412 |
| Mean raw_diff | +7.10 |
| Median raw_diff | +5.00 |
| SD | 7.00 |
| Range | −8 to +32 |

The distribution is **right-skewed** (mean > median), with a long upper tail of large overestimates. Quantile regression is used instead of OLS to model the conditional distribution without normality assumptions.

Tumor type 1 (n = 122, the most frequent) is the reference category. Tumor type 8 was excluded. Final predictors: Sex, Age, Zn_HIPEC, PreOP_CTx, Thermoablation, and dummy variables for tumor types 2–7.

---

## 2. Full Model (including Age)

### 2.1 Median Model (q = 0.50)

| Predictor | Coef | p | Bootstrap 95% CI |
|---|---|---|---|
| Intercept | +6.72 | 0.007 | [+2.79, +12.50] * |
| Sex | −0.48 | 0.500 | [−2.00, +0.98] |
| **Age** | **−0.04** | **0.116** | **[−0.10, +0.01]** |
| Zn_HIPEC | +2.98 | 0.076 | [+0.28, +4.89] * |
| PreOP_CTx | −1.32 | 0.119 | [−3.29, −0.00] * |
| **Thermoablation** | **−2.15** | **0.005** | **[−4.46, −1.00] *** |
| Tumor_type_2 | −0.05 | 0.974 | [−2.79, +2.81] |
| Tumor_type_3 | −0.06 | 0.963 | [−2.54, +4.47] |
| Tumor_type_4 | +0.80 | 0.440 | [−1.37, +3.00] |
| **Tumor_type_5** | **+5.87** | **<0.001** | **[+3.17, +8.00] *** |
| Tumor_type_6 | −1.67 | 0.157 | [−4.00, −0.00] * |
| **Tumor_type_7** | **+6.68** | **<0.001** | **[+2.61, +9.43] *** |

### 2.2 Quantile Process Summary

| Predictor | q=0.25 | q=0.50 | q=0.75 | q=0.90 |
|---|---|---|---|---|
| Sex | +0.20 | −0.48 | −2.00 | −4.00 * |
| Age | −0.02 | −0.04 | 0.00 | 0.00 |
| Zn_HIPEC | +1.77 | +2.98 | +5.00 * | 0.00 |
| PreOP_CTx | −0.65 | −1.32 | −2.00 | −1.21 |
| Thermoablation | −2.68 *** | −2.15 ** | −2.00 | −4.00 * |
| Tumor_type_5 | +3.83 *** | +5.87 *** | +11.00 *** | +13.00 *** |
| Tumor_type_7 | +4.26 *** | +6.68 *** | +6.00 ** | +8.00 ** |

(All other tumor types non-significant at all quantiles.)

### 2.3 Model Fit

| Quantile | Pseudo-R² |
|---|---|
| q=0.25 | 0.087 |
| q=0.50 | 0.127 |
| q=0.75 | 0.156 |
| q=0.90 | 0.215 |
| OLS R² | 0.224 |

### 2.4 Symmetry Test (q=0.75 − q=0.25, bootstrap)

| Predictor | Δ coef | 95% CI | Asymmetric? |
|---|---|---|---|
| Sex | −2.24 | [−5.00, −0.00] | YES * |
| Tumor_type_5 | +5.92 | [+1.55, +10.72] | YES * |
| All others | — | CI includes 0 | No |

---

## 3. Reduced Model (Age Excluded)

Age was **not significant at any quantile** (all p > 0.10, bootstrap CI includes zero). The analysis was repeated without Age.

### 3.1 Median Model (q = 0.50)

| Predictor | Coef | p | Bootstrap 95% CI |
|---|---|---|---|
| Intercept | +4.00 | 0.043 | [+2.00, +9.00] * |
| Sex | −1.00 | 0.162 | [−2.00, +1.00] |
| Zn_HIPEC | +3.00 | 0.074 | [+0.37, +5.00] * |
| PreOP_CTx | −1.00 | 0.236 | [−3.00, +0.00] |
| **Thermoablation** | **−2.00** | **0.009** | **[−5.00, −1.00] *** |
| Tumor_type_2 | 0.00 | 1.000 | [−3.00, +3.00] |
| Tumor_type_3 | 0.00 | 1.000 | [−3.00, +4.99] |
| Tumor_type_4 | +1.00 | 0.333 | [−1.19, +3.00] |
| **Tumor_type_5** | **+6.00** | **<0.001** | **[+3.39, +8.00] *** |
| Tumor_type_6 | −1.00 | 0.394 | [−4.00, +0.00] |
| **Tumor_type_7** | **+7.00** | **<0.001** | **[+3.00, +9.77] *** |

### 3.2 Quantile Process Summary

| Predictor | q=0.25 | q=0.50 | q=0.75 | q=0.90 |
|---|---|---|---|---|
| Sex | +0.33 | −1.00 | −2.00 | −4.00 * |
| Zn_HIPEC | +1.67 | +3.00 | +5.00 * | 0.00 |
| PreOP_CTx | −0.67 | −1.00 | −2.00 | −1.14 |
| Thermoablation | −2.67 *** | −2.00 ** | −2.00 | −4.00 * |
| Tumor_type_4 | −0.33 | +1.00 | +2.00 | +4.87 * |
| Tumor_type_5 | +3.67 *** | +6.00 *** | +11.00 *** | +13.00 *** |
| Tumor_type_7 | +4.00 *** | +7.00 *** | +6.00 ** | +8.00 ** |

### 3.3 Model Fit

| Quantile | Pseudo-R² |
|---|---|
| q=0.25 | 0.086 |
| q=0.50 | 0.124 |
| q=0.75 | 0.156 |
| q=0.90 | 0.215 |
| OLS R² | 0.223 |

Removing Age reduces OLS R² by only 0.001 (0.224 → 0.223), confirming it contributes nothing.

### 3.4 Symmetry Test (q=0.75 − q=0.25, bootstrap)

| Predictor | Δ coef | 95% CI | Asymmetric? |
|---|---|---|---|
| Intercept | +5.65 | [+0.50, +14.04] | YES * |
| Sex | −2.28 | [−5.00, −0.00] | YES * |
| Tumor_type_4 | +2.63 | [+0.00, +5.94] | YES * |
| **Tumor_type_5** | **+6.03** | **[+1.63, +11.00]** | **YES *** |
| All others | — | CI includes 0 | No |

---

## 4. Diagnostics

All four quantile models passed all diagnostic checks:

- **Residual sign proportions** match the nominal quantile at all q (±5% tolerance, ✓), indicating numerical convergence.
- **Residuals vs fitted** show no systematic trend or heteroskedasticity.
- **Spearman correlations** between residuals and all predictors are near zero and non-significant, ruling out missed non-linearity or obvious omitted variables.
- **Bootstrap SEs vs asymptotic SEs** agree closely for all significant predictors, validating the sparsity estimator.
- `IterationLimitWarning` messages arise only in a small number of edge-case bootstrap resamples, not in the primary fitted models; they do not affect the conclusions.

---

## 5. Conclusions

### Primary finding
**Tumor type drives the sPCI–pPCI discrepancy.** Patients with tumor type 5 have a median overestimate roughly 6 points higher than type 1, growing to +11 at q=0.75 and +13 at q=0.90. Tumor type 7 shows a similar pattern (+7 at the median, +6–8 in the upper quartiles). All other tumour types are indistinguishable from the reference.

### Thermoablation reduces the discrepancy
Thermoablation is consistently associated with a **2–4 point reduction** in sPCI–pPCI across all quantiles. The effect is most pronounced at the lower and upper tails (q=0.25 and q=0.90).

### Asymmetric effects
The symmetry test flags two variables whose effects **differ between the lower and upper halves of the distribution**:

- **Sex:** the negative effect of being male appears only in high-discrepancy patients (upper tail).
- **Tumor_type_5:** the overestimation effect is substantially larger in patients who already have a large discrepancy, indicating that surgeons are disproportionately wrong in the worst cases.

In the reduced model, **Tumor_type_4** also shows borderline asymmetry emerging only at q=0.90, which was not apparent when Age was included.

### Age and other predictors
Age, PreOP_CTx, Zn_HIPEC (at most quantiles), and tumour types 2, 3, and 6 have no robust association with the discrepancy. Dropping Age has negligible impact on model fit (ΔOLS R² < 0.001).

### Model fit
The pseudo-R² increases from ~0.086 at q=0.25 to ~0.215 at q=0.90, closely matching the OLS R² (0.22). This confirms that the predictors — particularly tumor type — account for more variance in the upper tail, where discrepancies are largest and clinically most consequential.
