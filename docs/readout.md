# Readout

All data are simulated, so every estimate can be checked against a known truth. The plan, with its rules, was committed before any estimator was run (`docs/design.md`); Part C was added as Amendment 1, again with its rules committed first. The code checks C1 (8 checks) and C2 passed (`results/check_log.txt`, `results/verdicts.txt`).

## Verdicts

| Rule | Result | Verdict |
|---|---|---|
| P1 | M3 (censored likelihood) bias: CL/F +2.0%, V/F +1.4% | met |
| P2 | Dropping BLQ (M1): V/F +7.8%; LLOQ/2 (M5): ka -5.9% | met |
| P3 | Gain from 100 to 200 mg, truth 7.1 points: one dose level 13.7, pooled 11.6, adjusted 7.2, by arm 7.2 | met |
| P4 | Chooses 200 mg when the gain is 7.1: pooled 61% vs adjusted 32%; when it is 17.1: adjusted 89% vs by arm 82% | met |
| P5 | Historical tumor model: arm effects 1.21 (truth 1.10) and -0.29 (-0.36), noise SD 0.121 (0.12), link 0.97 (0.80) | met |
| P6 | S1 (shrinks tumors, no survival benefit): response rule goes in 100% of studies in every cell, model rule in 4-8% | met |
| P7 | S2 (OS hazard ratio 0.71): model rule goes 27 to 44 points more often than the interim OS rule; 42 points at 60 per arm, 8-week scans | met |
| P8 | S2, 60 per arm: model rule goes in 86% with 6-week scans and 95% with 12-week scans | **not met** (opposite direction) |
| P9 | S4: median predicted hazard ratio 0.85 (truth 0.73); at 90 per arm the interim OS rule goes 10 to 16 points more often | met |
| C2 | Refitted M3 fits match the Part A objective values in 20 of 20 studies | met |
| P10 | 400 mg (not studied): median RMSE of the predicted adverse-event risk 9.1 points for the exposure-response model, 24.9 for XGBoost | met |
| P11 | 100 and 200 mg: 1.7 vs 6.4 and 2.7 vs 8.6 points | met |

## Part A: population PK and the dose choice

**Samples below the LLOQ.** 18.0% of the samples in the 20 studies were below the LLOQ, nearly all of them troughs. Means of the relative error over the 20 studies (SD of the estimates in brackets):

| Method | CL/F | V/F | ka | Albumin exponent |
|---|---|---|---|---|
| Censored, M3 | +2.0% (3.4) | +1.4% (2.6) | -0.9% (5.8) | +2.5% (24) |
| Set to LLOQ/2, M5 | +3.5% (3.6) | +0.6% (2.5) | -5.9% (5.5) | -12.0% (18) |
| Dropped, M1 | -5.4% (2.4) | +7.8% (2.8) | +10.5% (6.5) | -15.8% (25) |

Dropping the low troughs leaves the model to explain only the higher concentrations, so it stretches V/F and ka and lowers CL/F. LLOQ/2 did less harm here than dropping, but it still moved ka by 6%. The albumin exponent is the least precise estimate (SD about 20%), and both shortcuts pull it toward zero. None of the 60 fits failed; each took about a minute.

**Exposure-response and the dose choice.** Albumin lowers exposure (it raises clearance) and lowers the response rate, so patients with high exposure look better partly because they are healthier. In 2,000 studies with a true gain of 7.1 points from 100 to 200 mg:

| Estimator | Mean estimated gain | SD | Chooses 200 mg (gain at least 10) |
|---|---|---|---|
| Exposure-response, 200 mg arm only | 13.7 | 9.3 | 70% |
| Exposure-response, both doses | 11.6 | 5.8 | 61% |
| Arm comparison | 7.2 | 8.7 | 39% |
| Exposure-response, adjusted for albumin | 7.2 | 5.9 | 32% |

With a true gain of 17.1 points, the adjusted model chose 200 mg in 89% of studies and the arm comparison in 82%. The adjusted model is as unbiased as the randomized comparison and a third less variable, because it also uses the spread of exposures within each arm. The unadjusted models are about as precise but aim at the wrong number. Even the best estimator chooses 200 mg in about a third of studies when the gain is 7 points: 60 patients per arm cannot separate 7 from 10 points reliably.

The same pattern held with exposures from the empirical Bayes clearances of the M3 fits instead of the true ones (20 studies): pooled 12.1 and adjusted 7.2 points, against 11.4 and 7.0 with the true exposures.

## Part B: tumor dynamics and the go/no-go decision

**Historical model.** The Stein model fitted by SAEM to 600 patients in 26 seconds recovered the drug's effects on tumor shrinkage and regrowth and the measurement error (table in `results/historical_fit.csv`). The estimated link between regrowth rate and survival was 0.97 (SE 0.11) against a true 0.80, so the model-predicted hazard ratios lean toward benefit: in S2 the median prediction was 0.62 against a true 0.71.

**Decisions.** Each rule's threshold was set so that it goes in 10% of studies of a drug with no effect (S0). Go rates over all nine design cells (40, 60 or 90 per arm; scans every 6, 8 or 12 weeks):

| Scenario | True OS hazard ratio | Response rule | Model rule | Interim OS rule |
|---|---|---|---|---|
| S1 shrinks tumors only | 1.00 | 100% | 4-8% | 8-12% |
| S2 slows regrowth | 0.71 | 95-100% | 68-99% | 41-64% |
| S3 both | 0.71 | 100% | 71-99% | 37-61% |
| S4 part of the benefit outside the tumor | 0.73 | 100% | 16-43% | 35-53% |

- The response rate cannot tell S1 from S2: it goes almost always in both. The regrowth model goes in S1 no more often than for a drug with no effect, and in 92% of S2 studies at 60 per arm with 8-week scans.
- At week 50 most patients are still alive, so the interim OS hazard ratio is noisy. The model reads the same studies' tumor data and goes 27 to 44 points more often in S2.
- The model is only as good as its premise. In S4 about 60% of the hazard reduction (on the log scale) acts outside tumor growth; the model saw only the regrowth part (median prediction 0.85 against a true 0.73) and went less often than the interim OS rule.
- Fixed thresholds without calibration misbehave: a hazard ratio cut of 0.80 went in 16-29% of S0 studies for interim OS and 5-19% for the model. Thresholds should be set by simulation for each design.

**Scan interval (P8, not met).** The plan predicted that 6-week scans would beat 12-week scans. The opposite happened (86% vs 95% at 60 per arm in S2, and the same order in S3 and S4). Post hoc (`scripts/09_posthoc_scans.R`, 500 studies per cell at 60 per arm):

| Scans | Stop at progression | Patient scans after baseline | Follow-up per patient (weeks) | SD of the arm difference with no effect | S2 go |
|---|---|---|---|---|---|
| every 6 weeks | yes (as designed) | 2.3 | 14 | 0.22 | 84% |
| every 8 weeks | yes | 2.0 | 16 | 0.19 | 92% |
| every 12 weeks | yes | 1.5 | 17 | 0.17 | 95% |
| every 6 weeks | no | 4.5 | 27 | 0.13 | 99% |
| every 12 weeks | no | 2.0 | 24 | 0.14 | 99% |

Frequent scans call progression earlier, and scanning stops at progression, so each patient's tumor is followed for a shorter time. The regrowth rate is read from how the tumor grows over time, so shorter follow-up gives noisier estimates. When scanning continued after progression, the interval hardly mattered and power was near 99%. Tumor measurements taken after progression, where they exist, are worth keeping in tumor-dynamics models.

## Part C: adverse-event risk for a new patient

The true risk runs through each patient's exposure (AUC = dose / CL). Median RMSE of the predicted risk over 20 studies, in percentage points, against the true risk given the patient's weight and albumin:

| Dose | Exposure-response model | XGBoost | Logistic regression |
|---|---|---|---|
| 100 mg | 1.7 | 6.4 | 7.4 |
| 150 mg | 2.4 | 5.0 | 10.7 |
| 200 mg | 2.7 | 8.6 | 13.5 |
| 400 mg (not studied) | 9.1 | 24.9 | 21.9 |

At 400 mg the mean true risk was 37.4%. XGBoost predicted 14.2% on average: trees stay flat beyond the doses they were trained on, so they returned roughly the 200 mg risk. Logistic regression on log dose extrapolated the average well (38.2%) but spread its errors across patients, because with about 16 events per study it fitted noise in the five irrelevant covariates. The exposure-response model averaged 40.6% at 400 mg; with so few events its slope ranged from -0.1 to 3.2 across studies (truth 1.5), so its error at 400 mg is still sizeable (interquartile range 4 to 17 points).

This comparison favors the exposure-response model by construction: the truth has exactly its structure. What it shows is the value of that structure when it holds, above all for doses outside the studied range, not that machine learning cannot work with more data or better features.

## Limits

- Everything is simulated, and the fitted models have the true structure (one compartment, the Stein model, a logistic exposure-response curve). Real data add model misspecification, which this study does not test.
- A1 and Part C use 20 studies each, so their Monte Carlo SE is about 1 point for the PK biases.
- The tumor-survival link was estimated once, from one historical study; its error (0.97 against 0.80) carries into every predicted hazard ratio.
- Progression is RECIST-like but ignores new lesions and non-target disease; response is unconfirmed.
