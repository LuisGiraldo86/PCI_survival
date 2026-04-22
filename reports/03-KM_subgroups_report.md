# Kaplan-Meier Subgroup Survival Analysis

**Notebook:** `03-survival_KM_subgroups.ipynb`  
**Date:** April 21, 2026  
**Outcome:** Overall survival (OS), truncated at τ = 57 months  
**Dataset:** `GPT_processed_survival_data.csv`

---

## 1. Study Sample

| Parameter | Value |
|---|---|
| Total N | 264 |
| Death events | 101 (38.3%) |
| Truncation time τ | 57 months (inherited from `03-survival-kaplan_meier.ipynb`) |

Subgroup analyses were performed for **Sex**, **Tumor type**, and **Age** (continuous and tertile groups), all using truncated survival times and events.

---

## 2. Sex

### 2.1 Group Summary

| Sex | n | Events | Event rate |
|---|---|---|---|
| 0 (female) | 111 | 47 | 42.3% |
| 1 (male) | 153 | 54 | 35.3% |

Both groups meet the ≥ 10 events reliability threshold.

### 2.2 Log-rank Test

| Test statistic | p-value |
|---|---|
| 1.48 | 0.2234 |

No statistically significant difference in survival between sexes (p = 0.22).

### 2.3 Cox PH Model (univariable, Sex coded)

Proportional hazards assumption: **satisfied**.

The Cox model confirms no significant sex effect on survival. The HR estimate direction is consistent with the slightly higher event rate in females (Sex = 0), but the difference is not significant.

### 2.4 Power Assessment (Schoenfeld)

| Comparison | d | HR | Estimated power |
|---|---|---|---|
| Female vs. Male | 101 | 0.83 | 14.5% |

⚠ Power is well below the 80% threshold. The study is **underpowered to detect the observed sex effect** (HR ≈ 0.83). The non-significant log-rank result should be interpreted cautiously — a true moderate sex effect cannot be ruled out.

---

## 3. Completeness of Cytoreduction (CC)

### 3.1 Group Summary

| CC | n | Events | Event rate | Events >= 10? |
|---|---|---|---|---|
| 0 (complete) | 179 | 65 | 36.3% | OK |
| 1 (incomplete) | 84 | 35 | 41.7% | OK |

Both groups meet the >= 10 events reliability threshold. Total N = 263; total events = 100 (38.0%).

### 3.2 Log-rank Test

| Test statistic | p-value |
|---|---|
| 2.90 | 0.0888 |

No statistically significant difference in survival between CC groups at the 0.05 level (p = 0.089), though the trend is in the expected direction (higher event rate with incomplete cytoreduction).

### 3.3 Power Assessment (Schoenfeld)

| Comparison | d | HR | Estimated power |
|---|---|---|---|
| CC0 vs CC1 | 100 | 1.15 | 9.4% |

WARNING Power is well below the 80% threshold (9.4%). The study is **severely underpowered to detect the observed effect size** (HR ~ 1.15). The non-significant log-rank result is consistent with low power rather than a true null effect. Note that in multivariable analyses (see `04-cox_survival_report.md`), CC is a significant predictor of OS (HR = 1.67, p = 0.022) after adjusting for tumour type and age — the univariable signal is suppressed here by the small absolute HR and limited events per group.

---

## 4. Tumor Type

### 4.1 Group Summary

| Tumor | n | Events | Event rate | Events ≥ 10? |
|---|---|---|---|---|
| 1 | 99 | 48 | 48.5% | ✓ |
| 2 | 18 | 7 | 38.9% | ⚠ Low |
| 3 | 22 | 6 | 27.3% | ⚠ Low |
| 4 | 48 | 8 | 16.7% | ⚠ Low |
| 5 | 23 | 5 | 21.7% | ⚠ Low |
| 6 | 31 | 18 | 58.1% | ✓ |
| 7 | 13 | 4 | 30.8% | ⚠ Low |
| 8 | 10 | 5 | 50.0% | ⚠ Low |

Only tumor types 1 and 6 have ≥ 10 events; all other groups have sparse event counts.

### 4.2 Multivariate Log-rank Test

| Test statistic | p-value |
|---|---|
| 48.08 | 3.43 × 10⁻⁸ |

Tumor type is a **highly significant predictor of survival** (p < 0.001). The KM curves show clearly differentiated survival trajectories across tumor groups.

### 4.3 Pairwise Power Assessment (Schoenfeld)

Only **two pairwise comparisons** meet the ≥ 80% power threshold:

| Comparison | d | HR | Estimated power | Conclusion |
|---|---|---|---|---|
| T1 vs T4 | 56 | 0.34 | 97.9% | ✓ Sufficient power |
| T1 vs T5 | 53 | 0.45 | 83.1% | ✓ Sufficient power |
| T4 vs T6 | 26 | 3.48 | 88.9% | ✓ Sufficient power |

All other pairwise comparisons are underpowered (< 80%), primarily due to low event counts in most tumor-type subgroups. The overall multivariate test is robust (large combined events), but **individual pairwise contrasts should be interpreted with caution** except for the three adequately powered comparisons above.

The patterns are consistent with the Cox PH results from notebook `04-cox_standard.ipynb`: tumor types 3, 4, 5, and 7 show markedly lower event rates than type 1, while type 6 shows the highest event rate of all groups (58.1%).

---

## 5. Age

### 5.1 Tertile Groups (KM curves)

Age was divided into tertiles (Younger / Middle / Older) for visual inspection of KM curves. No formal test was performed on the tertile groups; the Cox PH model on continuous Age is the primary inferential analysis.

### 5.2 Cox PH Model (univariable, Age continuous)

| Parameter | Value |
|---|---|
| Total N | 264 |
| Total events | 101 (38.3%) |
| Age mean | 54.2 years |
| Age SD | 11.6 years |

| HR per 1-year increase | 95% CI | Wald p-value |
|---|---|---|
| 0.991 | [0.973, 1.008] | 0.287 |

Age is **not a significant predictor** of overall survival (p = 0.287). The HR per year is essentially 1.0.

Proportional hazards assumption: **satisfied**.

### 5.3 Power Assessment (Schoenfeld, per SD of Age)

| log(HR) per SD | Estimated power |
|---|---|
| 0.111 | 8.1% |

⚠ Power is extremely low (8.1%). The study has **insufficient power to detect the observed age effect** (HR ≈ 0.99 per year, equivalent to HR ≈ 0.90 per SD). A clinically meaningful age effect cannot be excluded on statistical grounds alone.

### 5.4 Predicted Survival Curves

Cox-predicted survival curves were plotted for the P25, P50, and P75 percentiles of Age, confirming the near-flat effect: the survival curves at different age values are very close across the follow-up period.

---

## 6. Summary

| Variable | Test | Test statistic | p-value | Significant | Power adequate? |
|---|---|---|---|---|---|
| Sex | Log-rank | 1.48 | 0.2234 | No | No (14.5%) |
| CC | Log-rank | 2.90 | 0.0888 | No (trend) | No (9.4%) |
| Tumor type | Multivariate log-rank | 48.08 | 3.43 × 10⁻⁸ | **Yes** | Overall yes; pairwise mostly no |
| Age | Cox Wald test | — | 0.287 | No | No (8.1%) |

**Key conclusion:** Tumor type is the only subgroup variable with a statistically significant association with overall survival in this univariable analysis. Sex, CC, and Age are all non-significant at the 0.05 level. However, all three are severely underpowered (power < 15%) for their observed effect sizes, so non-significance should not be interpreted as evidence of no effect. Notably, CC becomes a significant predictor in the multivariable Cox model (HR = 1.67, p = 0.022, see `04-cox_survival_report.md`), where tumour type is controlled for. The univariable CC log-rank result is therefore likely attenuated by confounding with tumour type distribution.
