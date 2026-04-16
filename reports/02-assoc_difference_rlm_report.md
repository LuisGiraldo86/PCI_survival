# Predictors of sPCI–pPCI Discordance: Robust Linear Model (RLM) Analysis

**Date:** April 16, 2026  
**Source:** `notebooks/02-assoc_difference_ols_robust.ipynb`

---

## 1. Background and Objective

Following the OLS analysis (see `02-assoc_difference_report.md`), which identified violations of normality and homoskedasticity, this analysis re-examines the predictors of sPCI–pPCI discordance using **Robust Linear Regression (RLM)** via iteratively reweighted least squares (IRLS) with a Huber-T M-estimator. This approach down-weights influential outliers and provides more reliable inference under non-normal, heavy-tailed residuals.

---

## 2. Study Sample

**N = 412** patients after excluding tumour type 8.  
Tumour type 1 (n = 122) serves as the reference category.

**Outcome:** `raw_diff` = sPCI − pPCI (signed discordance)

---

## 3. RLM Regression Results

**Model:** IRLS, Huber-T norm, scale estimated via MAD, covariance type H1  
**No. Observations:** 412 | **Df Residuals:** 400 | **Df Model:** 11 | **Iterations:** 26

### Regression Coefficients

| Predictor | Coefficient | SE (asymptotic) | SE (bootstrap) | z | p | 95% CI (bootstrap) |
|---|---|---|---|---|---|---|
| Intercept | +7.69 | 2.195 | 2.168 | 3.51 | **< 0.001** | [+3.66, +12.09] |
| Sex (male) | −0.83 | 0.635 | 0.654 | −1.30 | 0.193 | [−2.09, +0.41] |
| Age | −0.030 | 0.024 | 0.025 | −1.27 | 0.206 | [−0.079, +0.020] |
| Zn HIPEC | +2.54 | 1.487 | 1.206 | 1.71 | 0.088 | [−0.24, +4.45] |
| PreOP CTx | −1.41 | 0.750 | 0.682 | −1.89 | **0.044** | [−2.74, −0.02] |
| **Thermoablation** | **−2.75** | 0.682 | 0.821 | −4.03 | **< 0.001** | [−4.40, −1.27] |
| Tumour type 2 | −0.09 | 1.362 | 1.094 | −0.06 | 0.949 | [−2.16, +2.09] |
| Tumour type 3 | +0.65 | 1.208 | 1.183 | 0.54 | 0.593 | [−1.39, +3.23] |
| Tumour type 4 | +0.88 | 0.920 | 0.809 | 0.95 | 0.341 | [−0.61, +2.49] |
| **Tumour type 5** | **+5.74** | 0.821 | 0.990 | 6.99 | **< 0.001** | [+4.03, +7.92] |
| Tumour type 6 | −1.36 | 1.044 | 0.849 | −1.31 | 0.192 | [−2.94, +0.34] |
| **Tumour type 7** | **+5.26** | 1.205 | 1.502 | 4.36 | **< 0.001** | [+2.39, +8.25] |

**Weighted R²:** 0.246

### Key Findings

- **Tumour type 5** is the strongest predictor: +5.7 points more discordance than type 1 (p < 0.001), robust across all M-estimator norms tested.
- **Tumour type 7** shows +5.3 points more discordance than type 1 (p < 0.001), also highly consistent.
- **Thermoablation** is associated with −2.7 points less discordance (p < 0.001), the only continuous protective factor.
- **PreOP CTx** reaches borderline significance (p = 0.044, asymptotic; bootstrap CI excludes zero barely), associated with −1.4 points less discordance.
- Sex, Age, Zn HIPEC, and tumour types 2, 3, 4, and 6 are **not significant** predictors.

---

## 4. Model Fit and Residual Diagnostics

| Metric | Value |
|---|---|
| Weighted R² | 0.246 |
| Residual MAD | 3.36 |
| Residual RMSE | 6.22 |
| Outcome MAD | 4.00 |
| Outcome RMSE | 6.99 |
| Durbin-Watson | 2.014 (no autocorrelation) |
| Shapiro-Wilk (residuals) | W = 0.944, p < 0.001 (non-normal) |
| Skewness | 0.94 (right-skewed) |
| Kurtosis | 1.16 (leptokurtic) |

---

## 5. IRLS Weight Distribution

| Weight range | N | % |
|---|---|---|
| w = 1.0 (not downweighted) | 324 | 78.6% |
| 0.9 ≤ w < 1.0 | 14 | 3.4% |
| 0.5 ≤ w < 0.9 | 56 | 13.6% |
| 0.2 ≤ w < 0.5 | 18 | 4.4% |
| w < 0.2 | 0 | 0.0% |

The Huber threshold was 4.52 units; 153 residuals exceeded it. No observations received near-zero weights, indicating the sample contains influential but not extreme outliers.

---

## 6. Leverage and Influence

- **High-leverage observations** (h > 2p/n = 0.058): 26 (6.3%); max leverage = 0.111
- **Cook's D > 4/n**: 22 observations; **Cook's D > 1**: 0 observations; max Cook's D = 0.064

No observation exerts undue influence on the overall fit.

---

## 7. Sensitivity to M-estimator Norm

| Predictor | HuberT (default) | Bisquare (Tukey) | AndrewWave |
|---|---|---|---|
| PreOP CTx | −1.41 (p=0.044)* | −1.56 (p=0.027)* | −1.56 (p=0.026)* |
| Thermoablation | −2.75 (p<0.001)*** | −2.88 (p<0.001)*** | −2.88 (p<0.001)*** |
| Tumour type 5 | +5.74 (p<0.001)*** | +5.10 (p<0.001)*** | +5.07 (p<0.001)*** |
| Tumour type 7 | +5.26 (p<0.001)*** | +5.09 (p<0.001)*** | +5.09 (p<0.001)*** |
| Zn HIPEC | +2.54 (p=0.088) | +2.64 (p=0.047)* | +2.65 (p=0.049)* |
| Scale (MAD) | 5.097 | 4.994 | 4.990 |

All significant findings are **consistent across norms**. Zn HIPEC reaches marginal significance (p ≈ 0.047–0.049) under Bisquare and AndrewWave norms.

---

## 8. Standardised Residuals

| Threshold | N observations | % |
|---|---|---|
| \|z\| > 2 | 42 | 10.2% |
| \|z\| > 3 | 14 | 3.4% |
| \|z\| > 4 | 3 | 0.7% |
| Max \|z\| | 4.56 | — |

The proportion of large residuals is slightly elevated relative to a normal distribution, consistent with the right-skewed and leptokurtic outcome.

---

## 9. Comparison with OLS Results

| Finding | OLS | RLM | Conclusion |
|---|---|---|---|
| Tumour type 5 | +7.1*** | +5.7*** | Both significant; RLM attenuated (outlier influence reduced) |
| Tumour type 7 | +5.2*** | +5.3*** | Consistent |
| Thermoablation | −2.6*** | −2.7*** | Consistent |
| Sex (male) | −1.4* | −0.8 (NS) | OLS effect driven partly by influential observations |
| PreOP CTx | −0.9 (NS) | −1.4* | RLM reveals effect masked in OLS |
| R² / Weighted R² | 0.224 | 0.246 | Comparable overall fit |

The RLM analysis broadly confirms OLS findings but resolves the conflicting signal on Sex (no longer significant) and reveals a clearer effect for PreOP CTx.

---

## 10. Summary

The RLM model explains ~25% of the weighted variance in sPCI–pPCI discordance. Results are robust to outliers and consistent across multiple M-estimator norms. The dominant predictors are:

1. **Tumour type 5** (+5.7 points vs. type 1, p < 0.001) — highest intraoperative surprise
2. **Tumour type 7** (+5.3 points vs. type 1, p < 0.001) — second highest discordance
3. **Thermoablation** (−2.7 points, p < 0.001) — less discordance, likely reflecting patient selection
4. **PreOP CTx** (−1.4 points, p = 0.044) — modest protective association

Sex, age, Zn HIPEC, and tumour types 2, 3, 4, and 6 do not independently predict discordance. The residuals remain non-normal and right-skewed; quantile regression on the median is recommended as a further robustness check.
