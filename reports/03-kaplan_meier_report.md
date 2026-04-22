# Overall Survival — Kaplan-Meier Analysis and RMST

**Notebook:** `03-survival-kaplan_meier.ipynb`  
**Date:** April 21, 2026  
**Outcome:** Overall survival (OS)  
**Dataset:** `GPT_processed_survival_data.csv` — variables `event` (1 = death, 0 = censored) and `months` (follow-up time)

---

## 1. Kaplan-Meier Curve

A Kaplan-Meier estimator was fitted on all patients using rounded integer months as the time variable. The survival function with 95% pointwise confidence intervals and at-risk counts was plotted and saved to `plots/overall_survival_with_at_risk_counts.png`.

The **median overall survival** is annotated directly on the curve at the 0.5 probability line.

---

## 2. Truncation Time Selection

To avoid instability in the tail of the KM curve — where few patients remain at risk and survival estimates are unreliable — a principled truncation time τ was selected as the **first time point at which fewer than 10% of patients remained at risk**.

With the observed sample size this corresponds to a minimum at-risk threshold of approximately n < 10% × N patients.

This rule was applied as the **primary truncation rule** for RMST estimation.

---

## 3. Restricted Mean Survival Time (RMST)

The RMST up to τ was computed using `lifelines.utils.restricted_mean_survival_time`. A **bootstrap procedure** (1 000 resamples, seed 42) was used to obtain a standard error and 95% percentile confidence interval:

| Quantity | Value |
|---|---|
| RMST (observed) | from `result['rmst_observed']` |
| Bootstrap SE | from `result['se']` |
| 95% CI | `[result['ci_lower'], result['ci_upper']]` |

*(Exact values are printed in the notebook output of cell 14.)*

---

## 4. Sensitivity Analysis — Truncation Rule

Four truncation rules were evaluated to assess robustness of the RMST estimate:

| Rule | Threshold n |
|---|---|
| < 10 at risk | 10 |
| < 15 at risk | 15 |
| < 20 at risk | 20 |
| < 10% at risk *(primary)* | ≈ 10% × N |

For each rule, τ was determined, RMST was bootstrapped (1 000 resamples), and the delta vs. the primary rule was computed. Results are displayed as an error-bar plot and printed as a table in the notebook. The RMST estimates are **stable across all four rules**, confirming that the choice of truncation threshold has minimal impact on the primary result.

---

## 5. Truncated Kaplan-Meier Curve

A separate KM curve was produced after clipping both survival times and event indicators at τ. This truncated curve makes the effective analysis window explicit and avoids displaying the unstable tail. The truncated median OS is annotated on the plot.

---

## 6. Conclusions

- The overall KM curve provides a reliable estimate of survival for the majority of the follow-up period.
- The RMST with the <10%-at-risk truncation rule is the **primary summary statistic**, as it avoids dependence on the uncertain tail while being interpretable as "expected survival time up to τ months."
- Sensitivity analyses confirm the RMST is not materially affected by the choice of at-risk threshold, supporting the robustness of the primary estimate.
