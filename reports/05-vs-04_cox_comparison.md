# Cox PH Models: Comparison of Sex-Stratified vs. Sex-Excluded Analyses

**Notebooks compared:**
- `04-cox_standard.ipynb` — Sex used as stratification variable (separate baseline hazard per sex)
- `05-cox_no_sex.ipynb` — Sex dropped entirely (single shared baseline hazard)

**Dataset:** N = 263 patients, 100 events (deaths), follow-up truncated at 57 months.  
**Covariates (both):** Age, Tumor type (8 classes; T1 = reference), CC (complete cytoreduction), sPCI or pPCI.

---

## 1. Model Fit Summary

### Notebook 04 — Sex Stratified

| Model | C-index | Partial AIC | PCI HR (95% CI) | PCI p | LRT p |
|---|---|---|---|---|---|
| Base (Age + Tumor + CC) | 0.680 | 738.85 | — | — | — |
| Model A (+ sPCI) | 0.688 | 738.99 | 1.021 (0.991–1.052) | 0.171 | — |
| Model B (+ pPCI) | 0.698 | 736.58 | 1.046 (1.002–1.092) | **0.038** | — |

> Note: LRT p-values in notebook 04 appear inflated (≈0) due to a numerical issue in the comparison cell; p-values from the individual model summaries are used here.

### Notebook 05 — Sex Excluded

| Model | C-index | Partial AIC | PCI HR (95% CI) | PCI p | LRT p |
|---|---|---|---|---|---|
| Base (Age + Tumor + CC) | 0.688 | 880.66 | — | — | — |
| Model A (+ sPCI) | 0.697 | 880.49 | 1.023 (0.993–1.053) | 0.138 | 0.141 |
| Model B (+ pPCI) | 0.706 | 877.95 | 1.048 (1.005–1.094) | **0.030** | **0.030** |

> Note: AICs between the two notebooks are **not directly comparable** — stratification changes the partial log-likelihood scale (sex-specific baseline hazards absorb some residual variance).

---

## 2. Concordance (C-index) Comparison

| Model | C (Sex stratified) | C (Sex excluded) | Δ |
|---|---|---|---|
| Base | 0.680 | 0.688 | +0.008 |
| Model A (sPCI) | 0.688 | 0.697 | +0.009 |
| Model B (pPCI) | 0.698 | 0.706 | +0.008 |

Removing Sex from the model consistently raises the observed C-index by ~0.008 across all three models. This is counterintuitive but expected: when Sex is used as a stratum rather than a covariate (notebook 04), it does not contribute to the linear predictor and therefore cannot improve concordance on the original scale. In notebook 05, all residual variation previously absorbed by Sex stratification is left in the error term, paradoxically making the remaining predictors appear slightly sharper in global concordance — though absolute discrimination remains modest in both cases.

---

## 3. PCI Effect Estimates

| PCI type | HR (04, stratified) | HR (05, no Sex) | Δ HR |
|---|---|---|---|
| sPCI | 1.021 (0.991–1.052), p = 0.171 | 1.023 (0.993–1.053), p = 0.138 | +0.002 |
| pPCI | 1.046 (1.002–1.092), p = 0.038 | 1.048 (1.005–1.094), p = 0.030 | +0.002 |

The PCI effect estimates are **remarkably stable** across the two specifications — differences in HR are ≤ 0.002 and CIs are almost identical. This suggests Sex is not a meaningful confounder of the PCI–survival relationship.

- **sPCI** is non-significant in both analyses (p > 0.13).
- **pPCI** is borderline significant in both analyses, with every additional unit associated with ~4.6–4.8% higher hazard of death per unit increase.

---

## 4. Discrimination: Bootstrap ΔC-index (Model B − Model A)

| | Sex Stratified (04) | Sex Excluded (05) |
|---|---|---|
| Observed ΔC | +0.0099 | +0.0094 |
| Bootstrap 95% CI | [−0.0123, +0.0301] | [−0.0122, +0.0302] |
| Conclusion | No significant difference | No significant difference |

Bootstrap results are nearly identical. Neither analysis finds a statistically significant difference in discriminatory ability between sPCI and pPCI models. The confidence intervals straddle zero in both cases.

---

## 5. Other Covariates

Effects of tumour type and CC are consistent across both notebooks:

| Covariate | HR range (04) | HR range (05) | Notes |
|---|---|---|---|
| Age | ~0.99, NS | ~0.99, NS | Not significant in any model |
| CC | 1.40–1.67, p ≈ 0.02–0.18 | 1.38–1.65, p ≈ 0.02–0.19 | Significant in base model; attenuated when PCI added |
| T_3 | ~0.21–0.24, p < 0.002 | ~0.23–0.27, p < 0.002 | Strongly protective vs. T_1 |
| T_4 | ~0.14, p < 0.0005 | ~0.21, p < 0.0005 | Strongest protective effect (slightly attenuated without Sex) |
| T_5 | ~0.20–0.22, p < 0.002 | ~0.22–0.26, p < 0.003 | Protective vs. T_1 |
| T_6 | ~1.65–2.10, p ≈ 0.03–0.10 | ~1.68–2.20, p ≈ 0.03–0.08 | Increased hazard vs. T_1 |
| T_7 | ~0.28–0.34, p ≈ 0.02–0.04 | ~0.28–0.34, p ≈ 0.04–0.02 | Protective vs. T_1 |

The most notable difference is in **T_4**: the HR is ~0.14 with Sex stratification vs. ~0.21 without, suggesting Sex and T_4 share some variance. All other covariates are essentially unchanged.

---

## 6. Summary

1. **Sex has minimal impact on PCI effect estimates.** HRs for sPCI and pPCI differ by < 0.003 between the two model specifications.
2. **pPCI is the superior PCI measure** in both analyses: borderline significant (p ≈ 0.03–0.04), preferred by AIC (ΔAIC ≈ +2.4–2.5 over sPCI model), and marginally higher C-index — though the bootstrap ΔC is not statistically significant in either notebook.
3. **sPCI is not a significant predictor** of overall survival in either specification (p > 0.13).
4. **Tumour type dominates prognosis**, with T_4 carrying the most protective association vs. T_1.
5. Both analyses conclude that the difference in discrimination between sPCI and pPCI models is **not statistically significant** (bootstrap 95% CI crosses 0).
6. The sex-excluded specification (notebook 05) is preferred for comparability with interaction models (notebook 04-cox_interaction) that also exclude Sex.
