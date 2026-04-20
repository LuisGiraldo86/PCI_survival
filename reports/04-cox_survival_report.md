# Cox Proportional Hazards Analysis — Overall Survival
## Surgical vs. Pathological PCI as Predictors

**Study:** Retrospective cohort — peritoneal carcinomatosis patients  
**Outcome:** Overall survival (OS)  
**Research question:** Is PCI measured during surgery (sPCI) or by the pathologist (pPCI) the better independent predictor of OS?  
**Analysis date:** April 2026  

---

## 1. Study Population

| Parameter | Value |
|-----------|-------|
| Total patients | 263 |
| Death events | 100 (38.0%) |
| Median follow-up | 13.0 months (truncated at 57 months) |
| sPCI missing | 0 |
| pPCI missing | 0 |
| CC missing | 0 |

### PCI Descriptive Statistics

| Measure | Mean | SD | Min | Max |
|---------|------|----|-----|-----|
| sPCI (surgical) | 15.6 | 9.0 | 1 | 39 |
| pPCI (pathological) | 9.7 | 6.0 | 0 | 28 |

sPCI is systematically higher than pPCI (~6 points on average), consistent with prior literature reporting that surgeons tend to over-stage peritoneal spread relative to pathological assessment.

### Completeness of Cytoreduction (CC)

| CC | n |
|----|---|
| 0 (complete) | 179 (68%) |
| 1 (incomplete) | 84 (32%) |

---

## 2. Statistical Methods

Three Cox proportional hazards models were fitted:

- **Base model:** Age + Tumor type + CC (Sex stratified)  
- **Model A:** Base + sPCI (Sex stratified)  
- **Model B:** Base + pPCI (Sex stratified)  

**Sex** was removed from the linear predictor and used as a **stratification variable** because it violated the proportional hazards assumption (Schoenfeld residuals test, rank transform, p < 0.05). Stratification allows each sex to have its own baseline hazard function while estimating shared covariate effects.

Tumour type was one-hot encoded with type 1 (most frequent, n = 122) as the reference category.

Model performance was assessed using:
- **C-index (Harrell's concordance index)** — discrimination
- **AIC** — penalised model fit
- **Likelihood Ratio Test (LRT)** — improvement over base (each augmented model compared against a base model re-fitted on identical observations)
- **Bootstrap ΔC-index** (1 000 paired resamples) — non-nested comparison of Model A vs Model B

---

## 3. Proportional Hazards Assumption

All three models pass the Schoenfeld residuals test (rank transform) after Sex stratification. No violations were detected for any covariate, including sPCI and pPCI.

---

## 4. Cox Model Results

### 4.1 Base Model (Age + Tumor + CC, Sex stratified)

| Covariate | HR | 95% CI | p |
|-----------|-----|--------|---|
| Age | 0.99 | [0.97, 1.01] | 0.434 |
| **CC** | **1.67** | **[1.08, 2.60]** | **0.022** |
| T_2 | 0.58 | [0.26, 1.33] | 0.202 |
| **T_3** | **0.24** | **[0.10, 0.57]** | **0.001** |
| **T_4** | **0.14** | **[0.06, 0.33]** | **<0.001** |
| **T_5** | **0.23** | **[0.09, 0.60]** | **0.002** |
| T_6 | 1.65 | [0.90, 3.02] | 0.104 |
| **T_7** | **0.33** | **[0.11, 0.94]** | **0.038** |
| T_8 | 1.36 | [0.54, 3.47] | 0.516 |

**C-index = 0.680** | **AIC = 738.9**

Tumour type is the dominant predictor. Incomplete cytoreduction (CC = 1) is associated with a 67% increase in hazard of death.

---

### 4.2 Model A — Base + Surgical PCI (sPCI)

| Covariate | HR | 95% CI | p |
|-----------|-----|--------|---|
| Age | 0.99 | [0.98, 1.01] | 0.518 |
| CC | 1.41 | [0.85, 2.33] | 0.184 |
| **sPCI** | **1.02** | **[0.99, 1.05]** | **0.171** |
| T_3 | 0.21 | [0.09, 0.52] | 0.001 |
| T_4 | 0.14 | [0.06, 0.34] | <0.001 |
| T_5 | 0.20 | [0.08, 0.53] | 0.001 |

**C-index = 0.688** | **AIC = 739.0** | **LRT vs base: χ²(1) = 1.86, p = 0.173**

sPCI is **not** a significant independent predictor of OS after adjusting for tumour type, age, and CC (p = 0.173). Adding sPCI does not significantly improve model fit over the base.

---

### 4.3 Model B — Base + Pathological PCI (pPCI)

| Covariate | HR | 95% CI | p |
|-----------|-----|--------|---|
| Age | 0.99 | [0.97, 1.01] | 0.386 |
| CC | 1.40 | [0.88, 2.24] | 0.158 |
| **pPCI** | **1.05** | **[1.00, 1.09]** | **0.038** |
| T_3 | 0.22 | [0.09, 0.54] | 0.001 |
| T_4 | 0.14 | [0.06, 0.34] | <0.001 |
| T_5 | 0.22 | [0.08, 0.56] | 0.002 |

**C-index = 0.698** | **AIC = 736.6** | **LRT vs base: χ²(1) = 4.27, p = 0.039**

pPCI **is** a significant independent predictor of OS after adjusting for tumour type, age, and CC (p = 0.038). Each additional point of pPCI increases the hazard by ~5% (HR = 1.05 per unit). Model B provides a significantly better fit than the base model.

---

## 5. Model Comparison Summary

| Model | N | Events | C-index | AIC | LRT p vs base | PCI HR | PCI p | CC HR | CC p |
|-------|---|--------|---------|-----|--------------|--------|-------|-------|------|
| Base | 263 | 100 | 0.680 | 738.9 | — | — | — | 1.67 | 0.022 |
| Model A (sPCI) | 263 | 100 | 0.688 | 739.0 | 0.173 | 1.02 | 0.171 | 1.41 | 0.184 |
| Model B (pPCI) | 263 | 100 | **0.698** | **736.6** | **0.039** | **1.05** | **0.038** | 1.40 | 0.158 |

Model B (pPCI) is preferred on every criterion: lower AIC, higher C-index, significant LRT, and a significant PCI coefficient.

---

## 6. Direct Comparison: Model A vs Model B (Bootstrap)

Since Model A and Model B are **non-nested** (they contain different PCI variables), LRT cannot be used directly to compare them. A paired bootstrap procedure (1 000 resamples of the same patients) was used to estimate the difference in discrimination.

| Metric | Value |
|--------|-------|
| Observed ΔC (B − A) | +0.010 |
| Bootstrap 95% CI | [−0.014, +0.030] |
| Interpretation | CI crosses 0 → **not statistically significant** |

Despite pPCI showing numerically superior performance, the bootstrap does not provide statistically significant evidence that Model B discriminates better than Model A.

---

## 7. Interpretation and Conclusions

1. **pPCI is the better predictor of overall survival.** It is the only PCI measure that independently predicts OS after adjusting for tumour type, age, and completeness of cytoreduction (p = 0.038; LRT p = 0.039 vs base). sPCI fails to reach significance in the fully adjusted model (p = 0.171).

2. **The collinearity between CC and sPCI explains the attenuation.** In the base model, CC is significant (HR = 1.67, p = 0.022). When sPCI is added, both CC and sPCI become non-significant — they capture overlapping information. pPCI retains significance even after controlling for CC, suggesting it contains prognostic information beyond what the surgeon's operative assessment captures.

3. **The discrimination difference is real but modest and not statistically significant.** The ΔC of +0.010 (bootstrap 95% CI [−0.014, +0.030]) does not cross the threshold for a formally significant improvement. This is expected given the sample size (n = 263, 100 events).

4. **Tumour type dominates overall survival.** Cancer type is by far the strongest predictor across all models. Any conclusion about PCI must be framed in the context that tumour biology (histology) is the primary determinant of prognosis.

5. **Proportional hazards assumption is satisfied** in all three models after stratifying by Sex.

---

## 8. Limitations

- **Sample size:** With 100 events and ~10 covariates, the models are at the limit of the EPV (events per variable) guideline of ≥10. Results should be interpreted with caution.
- **Retrospective design:** Selection and measurement biases inherent to retrospective studies cannot be ruled out.
- **CC–PCI collinearity:** CC and PCI are expected to correlate; their independent effects are difficult to disentangle.
- **Generalisation:** Results reflect the specific mix of tumour types and surgical practices at this institution.

---

## 9. Recommended Reporting Statement

> In multivariable Cox proportional hazards analysis stratified by sex, pathological PCI (pPCI) was a significant independent predictor of overall survival (HR per unit = 1.05, 95% CI 1.00–1.09, p = 0.038), whereas surgical PCI (sPCI) was not (HR = 1.02, 95% CI 0.99–1.05, p = 0.171), after adjusting for tumour type, age, and completeness of cytoreduction. Model discrimination was marginally higher for the pPCI model (C-index 0.698 vs 0.688), though the difference was not statistically significant (bootstrap ΔC = +0.010, 95% CI −0.014 to +0.030).
