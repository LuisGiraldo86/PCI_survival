# Predictors of sPCI–pPCI Discordance: OLS Regression Analysis

**Date:** April 16, 2026  
**Source:** `notebooks/02-assoc_difference.ipynb`

---

## 1. Background and Objective

A systematic positive bias between intraoperative surgeon PCI (`sPCI`) and preoperative/pathological PCI (`pPCI`) has been established (mean sPCI − pPCI ≈ 6.9, see report `01-agreement_all_report.md`). This analysis investigates **which patient and clinical characteristics explain the magnitude of this discordance**, using ordinary least squares (OLS) regression.

---

## 2. Study Sample

After excluding tumour type 8, the analytic sample comprised **N = 412** patients.  
Tumour type 1 (n = 122) was set as the reference category.

**Predictors included:**

| Variable | Type | Notes |
|---|---|---|
| Sex | Binary | 0 = female, 1 = male |
| Age | Continuous | Years |
| Zn HIPEC | Binary | Zinc-based HIPEC regimen |
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
| SD | 7.00 | 6.68 |
| Median | 5.0 | 5.0 |
| Min | −8 | 0 |
| Max | 32 | 32 |

The positive mean and median of `raw_diff` confirm that surgeons consistently find **more** peritoneal disease intraoperatively than predicted preoperatively (median discordance of 5 PCI points). Both distributions depart significantly from normality (Shapiro-Wilk: p < 0.001 for both).

---

## 4. OLS Regression Results

**Outcome: `raw_diff` (sPCI − pPCI)**

| | F-statistic | df | p-value |
|---|---|---|---|
| Overall model | 10.52 | (11, 400) | 5.1 × 10⁻¹⁷ |

**Variance explained:** R² = 0.224, Adjusted R² = 0.203

The model is highly significant overall and explains approximately **20% of the variance** in discordance.

### Regression Coefficients

| Predictor | Coefficient | SE | t | p | 95% CI |
|---|---|---|---|---|---|
| Intercept | 7.23 | 2.44 | 2.96 | 0.003 | [2.44, 12.03] |
| **Sex (male)** | **−1.41** | 0.71 | −2.00 | **0.047** | [−2.79, −0.02] |
| Age | −0.02 | 0.03 | −0.69 | 0.492 | [−0.07, 0.03] |
| Zn HIPEC | 2.28 | 1.65 | 1.38 | 0.169 | [−0.97, 5.52] |
| PreOP CTx | −0.85 | 0.83 | −1.02 | 0.309 | [−2.49, 0.79] |
| **Thermoablation** | **−2.56** | 0.76 | −3.38 | **0.001** | [−4.05, −1.07] |
| Tumour type 2 | −0.08 | 1.51 | −0.06 | 0.955 | [−3.06, 2.89] |
| Tumour type 3 | 1.39 | 1.34 | 1.03 | 0.302 | [−1.25, 4.03] |
| Tumour type 4 | 1.44 | 1.02 | 1.41 | 0.159 | [−0.57, 3.45] |
| **Tumour type 5** | **+7.05** | 0.91 | 7.73 | **< 0.001** | [5.26, 8.85] |
| Tumour type 6 | −0.78 | 1.16 | −0.67 | 0.503 | [−3.06, 1.50] |
| **Tumour type 7** | **+5.24** | 1.34 | 3.92 | **< 0.001** | [2.61, 7.88] |

### Key Findings

- **Tumour type 5** is the strongest predictor: on average, these patients show **+7.1 points** more discordance than type 1 (p < 0.001).
- **Tumour type 7** shows **+5.2 points** more discordance than type 1 (p < 0.001).
- **Thermoablation** is associated with **2.6 points less** discordance (p = 0.001), possibly reflecting more localised and predictable disease in patients selected for this procedure.
- **Male sex** is associated with **1.4 points less** discordance than female sex (p = 0.047), a modest but statistically significant effect.
- Age, Zn HIPEC, PreOP CTx, and tumour types 2, 3, 4, and 6 are **not significant** predictors.

---

## 5. Residual Diagnostics

| Test | Statistic | p-value | Conclusion |
|---|---|---|---|
| Shapiro-Wilk (residuals) | W = 0.956 | < 0.001 | Residuals non-normal |
| Jarque-Bera | — | < 0.001 | Confirms non-normality |
| Breusch-Pagan (LM) | 39.88 | < 0.001 | Heteroskedasticity present |
| Breusch-Pagan (F) | 3.90 | < 0.001 | Heteroskedasticity present |
| Skewness | 0.83 | — | Right-skewed residuals |
| Kurtosis | 3.87 | — | Slightly leptokurtic |

Both OLS assumptions of **normality** and **homoskedasticity** of residuals are violated. The reported standard errors and p-values should be interpreted with caution; they are likely anti-conservative (i.e., p-values may be overstated).

---

## 6. Limitations and Recommendations

1. **Heteroskedasticity-robust standard errors** (HC3) should be computed to obtain reliable inference. This will not change the point estimates but may widen confidence intervals and inflate some p-values.
2. **Quantile regression** on the median would be more robust given the skewed, heteroskedastic outcome and would not rely on distributional assumptions.
3. The signed outcome (`raw_diff`) assumes direction matters. If only magnitude is of interest, rerunning the model on `abs_diff` should be considered.
4. Tumour types with small cell counts (types 2, 6, 7) warrant cautious interpretation due to limited statistical power.
5. The negative coefficient for **Thermoablation** likely reflects **selection bias** rather than a causal effect — patients undergoing thermoablation may have more resectable, well-demarcated disease that is easier to assess preoperatively.

---

## 7. Summary

The OLS model explains ~20% of variance in sPCI–pPCI discordance. The dominant drivers are **tumour type** (especially types 5 and 7, which carry far more intraoperative surprise than type 1), followed by **thermoablation** (associated with less discordance) and **sex** (males slightly less discordant). Age and preoperative chemotherapy do not independently predict discordance. OLS assumptions are violated and robust inference methods are advised before drawing definitive conclusions.
