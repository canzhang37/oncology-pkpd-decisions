# Plan: population PK/PD and tumor-dynamics decisions for an oncology dose

Written before any estimator was run on the simulated studies. The truths in `results/truth_*.csv` were computed first, only to set the scenario parameters. Results are read only against the rules below; anything added later is labeled as an amendment.

## Questions

A fictional oral kinase inhibitor is in a randomized phase 2 study. Two kinds of decisions follow from such a study, and both rest on models fitted to sparse, imperfect data.

**Part A, dose.** Patients get 100 mg or 200 mg once daily.

- A1. About one sample in six is below the assay's lower limit of quantification (LLOQ). How much do the usual shortcuts, dropping those samples (Beal's M1) or setting them to LLOQ/2 (M5), bias a population PK model, compared with treating them as censored (M3)?
- A2. Albumin affects both clearance and the chance of response. How far does an exposure-response analysis overstate the gain from the higher dose, and does that change the choice between 100 and 200 mg?

**Part B, go/no-go.** The drug is compared with control in a randomized phase 2 study read 50 weeks after the first patient starts.

- B1. Can a decision based on the modeled tumor regrowth rate tell a drug that only shrinks tumors from one that improves survival, when the response rate cannot?
- B2. How do the number of patients and the scan interval change these decisions?

## Part A: drug, data and methods

PK model (simulation truth): one compartment with first-order absorption, 100 or 200 mg every 24 hours.

| Parameter | Typical value | Between-patient SD (log scale) |
|---|---|---|
| Absorption rate ka | 1.0 per hour | 0.4 |
| Clearance CL/F | 10 L/h × (WT/70)^0.75 × (ALB/40)^-1.0 | 0.3 |
| Volume V/F | 150 L × (WT/70) | 0.25 |

Weight is log-normal around 70 kg (SD 0.2 on the log scale, limited to 40-130 kg); albumin is normal with mean 38 and SD 5 g/L, limited to 25-50. Residual error is proportional 15% plus additive 0.01 mg/L. LLOQ is 0.2 mg/L. Samples are taken 0.5, 2 and 6 hours after the first dose, just before and 2 hours after the day 15 dose, and just before the day 29 dose. In large simulated populations 17.5% of samples are below the LLOQ, mostly the troughs (`results/truth_blq.csv`).

Each study has 60 patients per dose. The data sets for nlmixr2 write the doses before each trough as steady-state records; by day 15 the drug is within 0.5% of steady state for every simulated patient (check C1a2), so this changes the model by much less than the residual error.

**A1 fits.** The fitted model has the true structure (weight effects fixed at their allometric values, albumin exponent estimated) and is fitted by FOCEi in nlmixr2 under three treatments of BLQ samples: M3 (censored likelihood), M1 (dropped) and M5 (set to LLOQ/2). Start values are off from the truth by 20-50%. No covariance step is run (about 50 seconds per fit without it). R = 20 replicate studies (60 fits). For each method: mean relative bias of CL/F, V/F, ka and the albumin exponent, the SD of the estimates, and the share of fits that fail.

**A2 exposure-response.** Objective response follows

logit P(response) = b0 + b_exp × log(Cavg / Cref) + 0.6 × (ALB - 38) / 5,

where Cavg = dose / (24 × CL/F) is the average steady-state concentration, Cref the typical Cavg at 100 mg, and b0 gives a 30% response rate at Cref. Low albumin raises clearance (lower exposure) and lowers the chance of response, so exposure is confounded with prognosis. Two truths:

| Truth | b_exp | Response at 100 mg | Gain from 100 to 200 mg |
|---|---|---|---|
| small | 0.5 | 31.1% | 7.1 points |
| large | 1.2 | 31.1% | 17.1 points |

(`results/truth_er.csv`; the gain is the population average from 1,000,000 simulated patients.)

The estimand is the population gain in response rate from 100 to 200 mg. Estimators, each applied to a study of 60 patients per dose:

1. **One dose level**: logistic regression of response on log Cavg in the 200 mg arm only, as when a study tests one dose; the gain is predicted for those patients at 100 and 200 mg (g-computation over their clearances).
2. **Pooled**: the same model on both arms.
3. **Adjusted**: pooled, with albumin in the model.
4. **By arm**: the difference in observed response rates between the arms.

The dose rule: choose 200 mg only if the estimated gain is at least 10 points, otherwise 100 mg. 2,000 studies per truth, with exposures known exactly. As a realistic check, the 20 A1 studies also get responses from the small truth, and estimators 2 and 3 are refitted with Cavg from the empirical Bayes clearances of the M3 and M1 fits (descriptive only).

## Part B: tumor, survival and the decision

Tumor size (sum of longest diameters) follows the Stein model, TS(t) = TS0 × (exp(-d t) + exp(g t) - 1), with t in weeks. Control: TS0 log-normal around 70 mm (SD 0.5, limited to 20-250 mm), shrinkage rate d around 0.02 and regrowth rate g around 0.025 per week, each with SD 0.6 on the log scale. Measurements carry log-normal error with SD 0.12. Scans run every 6, 8 or 12 weeks from randomization until progression (20% and 5 mm above the smallest size so far), death or the data cutoff. A responder has a scan at least 30% below baseline before progression (unconfirmed).

Overall survival is Weibull (shape 1.3, typical control median 52 weeks) with log hazard 0.8 × log(g / 0.025) + 0.3 × log(TS0 / 70), plus, in one scenario, a drug effect outside the tumor.

| Scenario | Drug effect | Response rate, control / drug | True OS hazard ratio |
|---|---|---|---|
| S0 | none | 7.6% / 7.7% | 0.99 |
| S1 | shrinks tumors (d × 5), regrowth unchanged | 8.1% / 55.3% | 1.00 |
| S2 | slows regrowth (g × 0.6, d × 1.5) | 7.5% / 31.8% | 0.71 |
| S3 | both (d × 5, g × 0.6) | 7.8% / 75.4% | 0.71 |
| S4 | d × 5, g × 0.85, and hazard × 0.8 outside the tumor | 8.1% / 62.8% | 0.73 |

(`results/truth_tgi.csv`; response rates with 8-week scans over two years, hazard ratios from Cox fits to 20,000 uncensored patients per arm.)

**Historical study (B0).** The model behind the decision is built on an earlier randomized study: 300 patients per arm, control vs. an earlier drug (d × 3, g × 0.7), scans every 8 weeks, cutoff 130 weeks. The Stein model with arm effects on log d and log g and random effects on TS0, d and g is fitted by SAEM in nlmixr2. Then each patient's posterior mean of log g is computed on a grid (40 × 40 values of log d and log g; baseline size profiled out) using the fitted spread and measurement error, and a Cox model of OS on that and log baseline size gives the link coefficient b_g.

**Phase 2 studies.** Randomized 1:1, 40, 60 or 90 patients per arm, accrual uniform over 26 weeks, data cutoff at week 50 (the last patient has at least 24 weeks). 1,000 studies per scenario and design cell (5 × 3 × 3 cells, 45,000 studies). Each study gives three metrics:

1. **Response rate difference**, drug minus control.
2. **Model-predicted hazard ratio** exp(b_g × (mu_g,drug - mu_g,control)), where each arm's mean log regrowth rate mu_g is estimated by marginal maximum likelihood on the grid with the spread and noise fixed at the historical fit. Patients without a scan after baseline do not contribute.
3. **Interim OS hazard ratio** from a Cox model at the cutoff.

Decision rules:

- **Fixed thresholds** (descriptive): go if the response difference is at least 15 points; if the predicted hazard ratio is at most 0.80; if the interim OS hazard ratio is at most 0.80.
- **Calibrated thresholds** (the comparison the rules below use): in each design cell, each metric's threshold is set so that 10% of the S0 studies go. Go rates in S1 to S4 are then compared at the same 10% false-go rate under no effect.

## Rules, fixed before the results

C1 (code checks in `scripts/01_check.R`, must pass before any result is read):
- (a) the closed-form concentrations equal rxode2 solutions with all 29 doses within 1e-6 (relative); (a2) the steady-state records in the fitted data sets give concentrations within 0.5% of the full dosing history for 1,000 patients;
- (b) on noise-free tumor curves the grid's best point is within one grid step of the true log d and log g;
- (c) for 3,000 control patients with 6-week scans to week 60, the estimated mean log g is within 0.05 of the truth;
- (d) a Cox fit to 50,000 simulated survival times recovers 0.8 and 0.3 within 0.03;
- (e) M3 data sets keep every BLQ sample as censored at the LLOQ and M1 data sets drop them;
- (f) in a study of 50,000 patients per arm, the adjusted estimator is within 0.5 points and the by-arm estimator within 1 point of the true gain.

A1 (20 studies; the Monte Carlo SE of a mean relative bias is about 1-2%):
- P1. M3: mean relative bias within ±5% for CL/F and V/F.
- P2. M1 (dropping BLQ) biases CL/F or V/F by more than 5%; M5 (LLOQ/2) biases at least one of CL/F, V/F and ka by more than 5%.

A2 (2,000 studies per truth; Monte Carlo SE of a mean gain about 0.2 points, of a choice rate about 1 point):
- P3. Small truth: the pooled estimator without albumin overstates the gain by more than 30% (mean above 9.2 points), and the one-dose-level estimator overstates it more than the pooled one. The adjusted and by-arm estimators are within 1.5 points of the truth on average.
- P4. Small truth (the right choice is 100 mg): the pooled estimator without albumin chooses 200 mg at least 1.5 times as often as the adjusted one. Large truth (the right choice is 200 mg): the adjusted estimator chooses 200 mg more often than the by-arm comparison.

B0:
- P5. The SAEM fit recovers the arm effects on log d and log g within 0.15 and the measurement SD within 20%; the estimated link b_g is within 0.25 of 0.8.

B (calibrated thresholds; Monte Carlo SE of a go rate is at most 1.6 points):
- P6. S1 (shrinks tumors, no survival benefit): the response rule goes in more than 80% of studies in every design cell; the model rule in less than 25%.
- P7. S2 (hazard ratio 0.71 through slower regrowth): the model rule goes more often than the interim OS rule in every design cell, and by at least 10 points at 60 patients per arm with 8-week scans.
- P8. S2, 60 patients per arm: the model rule's go rate is at least 5 points lower with 12-week scans than with 6-week scans.
- P9. S4 (much of the benefit outside tumor regrowth): the model rule's median predicted hazard ratio is above 0.80, and with 90 patients per arm the interim OS rule goes more often than the model rule at every scan interval.

The fixed-threshold rules and S3 are reported without a rule.

## Amendment 1: adverse-event risk for a new patient (Part C)

Added after Parts A and B were planned and before any Part C code was run. Parts A and B are unchanged.

**Question.** A clinical pharmacology team is asked for a new patient's risk of a grade 3+ adverse event at a given dose, including a dose above the range studied. Does an exposure-response model built on the population PK fit do better than a flexible machine learning model trained on dose and covariates?

**Truth.** Each patient's risk depends only on their own steady-state exposure: P(grade 3+ AE) = logistic(-1.5 + 1.5 log(AUC / 20)), with AUC = dose / CL (mg·h/L). Dose, weight and albumin act only through clearance. For a typical patient the risk is 7.3% at 100 mg, 18.2% at 200 mg and 38.7% at 400 mg.

**Training data.** The 20 A1 studies (120 patients each, 100 or 200 mg), rebuilt from the same seeds; each patient also gets five covariates with no role in the truth (age, sex, ECOG 0/1, creatinine clearance, LDH) and an adverse-event outcome. The M3 fit of each study is refitted to recover its empirical Bayes clearances.

**Methods.**
1. **Exposure-response model**: logistic regression of the adverse event on log AUC, with each training patient's AUC from their empirical Bayes clearance. A new patient has no PK samples, so their risk is the fitted curve averaged over the fitted distribution of clearance given their weight and albumin (20-point Gauss-Hermite quadrature).
2. **XGBoost**: binary logistic on dose, weight, albumin and the five noise covariates; depth 1, 2 or 3 and the number of rounds (learning rate 0.05, early stopping after 50 rounds) chosen by 5-fold cross-validated log loss.
3. **Logistic regression** on log dose, log weight, log albumin and the noise covariates.

**Test set.** 5,000 new patients at each of 100 and 200 mg (studied), 150 mg (between) and 400 mg (above the studied range), per study. Metric: root mean squared error against the oracle risk, the true curve averaged over the true distribution of clearance given weight and albumin.

**Rules.**
- C2. The refitted M3 fits reproduce the A1 objective function values within 0.001.
- P10. At 400 mg, the exposure-response model's median RMSE over the 20 studies is less than half of XGBoost's.
- P11. At 100 and 200 mg, the exposure-response model's median RMSE is no larger than XGBoost's.

Logistic regression and the 150 mg results are reported without a rule.
