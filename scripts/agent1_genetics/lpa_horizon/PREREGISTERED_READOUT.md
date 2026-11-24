# Genetic emulation of Lp(a)HORIZON: pre-registered readout

Written 2026-11-20, before any estimate relating Lp(a) genotype to a later event in
people with established disease exists in this project. The code that implements these
rules is `lpa_horizon.R` in this folder (`verdict()` is the section "Verdict" below).

## Question

Genetically higher Lp(a) raises the risk of a first vascular event. Lp(a)HORIZON lowered
Lp(a) by a time-averaged 144 nmol/L for a median 4.7 years in people with established
disease and found no benefit (hazard ratio 1.02, 95% CI 0.90 to 1.15; Cho et al., NEJM
2026, doi 10.1056/NEJMoa2610159). Does the genetic effect of Lp(a) persist once disease
is established, which is the population the trial treated?

## What had been seen when this was written

- Counts and exposure distributions only (`feas_counts.R`, LS6 job 3600440): 16,955
  participants with established disease and 4,790 later events; predicted-Lp(a) SD
  59.6 nmol/L. That job never tabulated an outcome against Lp(a) or genotype.
- The published result for the paper's Lp(a) score on incident CAD in people free of CAD
  (hazard ratio 1.096 per SD, 1.083 to 1.110; 379,896 participants, 32,053 events).
- Published work by others: Zewinger et al. 2017 (doi 10.1016/S2213-8587(17)30096-7), no
  association of two *LPA* variants with mortality in 106,353 patients with established
  coronary disease; Patel et al. 2019 (doi 10.1161/CIRCGEN.119.002471), 9p21 odds ratio
  1.02 for subsequent coronary death or infarction against 1.20 for first events.

## Data and definitions

Cohort: the manuscript's unrelated European-ancestry set with *LPA* dosages and a
baseline assessment date (402,119). Follow-up ends at death, loss to follow-up or
2022-10-31.

Hospital data on LS6 hold the first date of each distinct code only. A later event of
the same type is therefore visible only under a different 4-character code; new event
types, new procedure codes and deaths are complete.

| Cohort | Definition |
|---|---|
| A, established disease (primary) | Before baseline: myocardial infarction (ICD-10 I21 to I23, I25.2; first-occurrence and algorithmic fields; self-report), ischaemic stroke (I63, I64; same sources) or peripheral artery disease (I70.2, I73.8, I73.9, I74.3 to I74.5; self-report). These are the trial's three qualifying conditions. |
| C, primary prevention | None of A, nor angina, other chronic coronary disease, coronary revascularization, TIA or unspecified stroke before baseline. |
| C standardized | C weighted to the age band (under 45, 45 to 49, 50 to 54, 55 to 59, 60 to 64, 65 and over) by sex mix of A. This is the comparator for the primary readout. |
| D, replication | Members of C with a first infarction or ischaemic stroke during follow-up who survive 30 days, followed from day 30. Disjoint from A. |
| B, sensitivity | Any vascular disease before baseline (A plus the conditions excluded from C). |

Endpoint: first of cardiovascular death (underlying cause ICD-10 chapter I), myocardial
infarction (I21, I22), ischaemic stroke (I63, I64) or coronary revascularization (OPCS-4
K40 to K46, K49, K50, K75) after time zero. Urgent and elective revascularization cannot
be told apart. Primary analysis uses all follow-up; a five-year horizon is the
trial-matched secondary analysis.

Exposures, all oriented so that a positive coefficient means higher risk:

- **Lp(a):** genotype-predicted Lp(a) in nmol/L from a censored-normal model (left-censored
  at 3.8, right-censored at 189, using the reportability field) fitted jointly on the 45
  *LPA* cis variants. Coefficients are per 144 nmol/L higher; results are also printed per
  144 nmol/L lower, the trial's direction.
- **LDL-C (positive control):** the paper's 167-variant LDL-C score, per SD. Lowering LDL-C
  after a first event is proven to reduce later events, so a design that cannot see this
  score in cohort A cannot be read for Lp(a).
- **Non-lipid CAD polygenic score (comparator, not in the verdict):** the UK Biobank-free
  CAD score restricted to variants more than 500 kb from any ApoB locus, per SD.

Model: Cox regression on time since time zero with the three exposures together, age,
sex and ten ancestry principal components. Cohort C standardized uses the weights above
with a robust variance.

## Gates (no estimate is read if any fails)

1. Cohort A size, its event count and cohort C size equal the feasibility job's counts.
2. Exposure model: no variant dropped, rs10455872 between 100 and 120 nmol/L per allele,
   r2 above 0.45 against in-range measurements.
3. The published Lp(a) score model reproduces: 379,896 participants, 32,053 events,
   hazard ratio 1.096 (1.083 to 1.110) to three decimals.
4. The synthetic self-test passes in full (data-assembly fixture, censored model,
   planted effects, every verdict branch).

Status of gate 4 at the first commit: LS6 job 3600625, 46 checks passed, script md5
`6e1c268df17113e7959d4078ee830b64`. Planted against recovered log hazard ratios for Lp(a)
in established disease: 0.00 against 0.042 (SE 0.039), 0.25 against 0.193 (SE 0.037),
0.00 against -0.007 (SE 0.039); the three scenarios returned the three intended verdicts.

## Primary readout

For each exposure, b_A is the log hazard ratio in cohort A and b_C in cohort C
standardized. A and C share no participants, so their estimates are independent.

- **Anchor:** b_C above 0 with one-sided P below 0.001, for Lp(a) and for LDL-C.
- **Present:** b_A above 0 with one-sided P below 0.025.
- **Attenuated:** b_A below b_C with two-sided P below 0.05 for the difference.
- **Specific:** the ratio b_A / b_C is smaller for Lp(a) than for LDL-C, one-sided P
  below 0.05 (delta method).

## Verdict

| Lp(a) in established disease | LDL-C control present in A | Lp(a) more attenuated than LDL-C | Verdict |
|---|---|---|---|
| anchor fails for either exposure | | | NO_ANCHOR |
| present, not attenuated | any | any | GENETICS_DOES_NOT_REPRODUCE_TRIAL |
| present, attenuated | any | any | PARTIAL_PERSISTENCE |
| not present, not attenuated | any | any | INDETERMINATE |
| not present, attenuated | no | any | UNINFORMATIVE_DESIGN |
| not present, attenuated | yes | yes | STAGE_SPECIFIC_FOR_LPA |
| not present, attenuated | yes | no | ATTENUATED_NOT_SPECIFIC |

How each is to be written up:

- **STAGE_SPECIFIC_FOR_LPA:** Lp(a) genotype predicts first events but not later ones,
  while LDL-C genotype predicts both. Genetics agrees with the trial when asked the
  trial's question.
- **ATTENUATED_NOT_SPECIFIC:** the Lp(a) effect is gone and the LDL-C effect is still
  detectable, but the data do not show that Lp(a) lost more. Report as consistent with the
  trial, not as specific to Lp(a).
- **UNINFORMATIVE_DESIGN:** a proven-treatable exposure is also invisible, so selection
  on having survived a first event, or loss of first-event genetics in general, explains
  the result. Nothing is concluded about Lp(a).
- **PARTIAL_PERSISTENCE** and **GENETICS_DOES_NOT_REPRODUCE_TRIAL:** the genetic effect
  persists, so the trial result is not explained by disease stage; treatment duration or
  the intervention remain.
- **INDETERMINATE** and **NO_ANCHOR:** no conclusion.

A flag `bias_material` is added to any verdict when the selection bias implied by
measured factors left out of the primary model (BMI, smoking, diabetes, hypertension)
exceeds a quarter of b_C for Lp(a). The bias is the sum over factors of (slope of the
factor on predicted Lp(a) in A minus the same slope in C) times the factor's log hazard
ratio in A. Unmeasured factors cannot be bounded this way and this is a stated limit.

## Secondary analyses (reported in full, no verdict)

Five-year horizon in cohort A, as a hazard ratio per 144 nmol/L lower with its
confidence interval and a test against the trial's estimate; cohorts D and B; A and D
pooled; cohort C unweighted; prior infarction, stroke and peripheral disease separately;
three-point composite and each component; one exposure per model; the paper's Lp(a) score
in place of predicted Lp(a); adjustment for the four risk factors; statin users in A;
proportional-hazards test for Lp(a) in A.

## Limits fixed in advance

- Genotype differences are lifelong, including before the first event. A null in cohort A
  is therefore stronger evidence than the trial's; a positive result would not separate
  treatment duration from the intervention.
- Most of cohort A has Lp(a) below the trial's entry range (1,845 have a measured value
  above 189 nmol/L), so the estimate is a per-nmol/L effect across the whole range. A
  threshold effect at very high Lp(a) cannot be excluded.
- Later infarctions are under-counted (first date per code), cohort A was less
  intensively treated than the trial's patients (median LDL-C 2.7 against about 1.7
  mmol/L), and participants are of European ancestry.
- No external replication cohort is available.

## Prediction, written before the run

Lp(a) will be attenuated in cohort A (the 2017 mortality result points that way). I am
unsure the LDL-C control will be detected, because nearly four in five of cohort A take
cholesterol-lowering drugs. The most likely verdicts are ATTENUATED_NOT_SPECIFIC or
UNINFORMATIVE_DESIGN; STAGE_SPECIFIC_FOR_LPA would be the most useful and I give it
about one chance in four.

## Deviations

Any change after this file is committed is listed here with its date, its reason and
whether an estimate had been seen.

- None at the time of the first commit.

## Outcome of the pre-registered run

LS6 job 3600628 on 2026-11-20, code at commit 238b4c3 (script md5
`6e1c268df17113e7959d4078ee830b64`). All 13 gates passed. Verdict of record:
**PARTIAL_PERSISTENCE**. Lp(a) per 144 nmol/L higher: cohort A hazard ratio 1.19 (1.12 to
1.27), cohort C standardized 1.35 (1.30 to 1.40), difference P 0.00094, ratio of log
hazard ratios 0.58 (0.35 to 0.81). LDL-C control present in A (1.07, 1.04 to 1.10). Lp(a)
not more attenuated than LDL-C (P 0.24). Bias flag not raised. The written prediction
(attenuated to absent) was wrong. Outputs: `exports/agent1/lpa_horizon/emulation/`.

## Post hoc analyses (2026-11-20, after the outcome above had been read)

Nothing here changes the verdict of record. Both analyses are in
`lpa_horizon_posthoc.R`, which sources the pre-registered script without editing it, and
both were written after seeing the estimates, so they are exploratory.

One piece of arithmetic was already done by hand after the run and is post hoc too: the
three-point composite at five years (0.88 per 144 nmol/L lower, 0.78 to 1.00) is not
significantly different from the trial's 1.02 (P about 0.10).

**P1. Qualifying condition.** Seen in the results: cohort A pools prior infarction, prior
ischaemic stroke and peripheral artery disease, the primary model does not adjust for
which, and the Lp(a) association with later stroke ran the wrong way (1.19 per 144 nmol/L
lower, 1.04 to 1.37). If high Lp(a) brings people into the cohort through infarction
rather than stroke, a pooled model can show this without any real effect. Analysis: (a)
mean predicted Lp(a) by qualifying condition (no outcome); (b) the primary model in A
stratified on the seven combinations of the three conditions, every endpoint and both
horizons, with the stage contrast recomputed; (c) every endpoint within prior infarction,
prior stroke, peripheral disease, infarction without stroke and stroke without
infarction.

**P2. What trial effect does the genetic estimate predict?** Genotype differences are
lifelong and a trial lasts about five years, so a genetic hazard ratio cannot be compared
with a trial's directly. LDL-C supplies the conversion because both are known for it.
Analysis: in each cohort, (a) the difference in measured LDL-C per SD of the LDL-C score,
adjusted as in the outcome model; (b) the LDL-C difference with the same genetic effect
on the four-point composite as 144 nmol/L of Lp(a), with a Fieller interval; (c) the
predicted hazard ratio for a trial lowering Lp(a) by 144 nmol/L, as 0.78 raised to that
equivalence. The constant 0.78 (0.76 to 0.80) per 1 mmol/L is the Cholesterol Treatment
Trialists' rate ratio for major vascular events (coronary death, infarction, coronary
revascularization, ischaemic stroke) over a median of about five years (Lancet 2010, doi
10.1016/S0140-6736(10)61350-5); its uncertainty is ignored. Cohort C is run as a check on
the method: it should give roughly the published equivalence of about 0.6 mmol/L per 144
nmol/L (Burgess et al. 2018). Assumption: per unit of genetic effect, short-term benefit
is the same fraction of lifelong benefit for Lp(a) as for LDL-C. Stated weakness: the
LDL-C difference between genotypes is smaller on statins than it was before them, so
measured LDL-C in cohort A understates the lifelong difference; the calibration is
repeated with LDL-C divided by 0.7 in statin users to show the sensitivity.

## Outcome of the post hoc analyses

Run 2026-11-20 at commit `83f27b3` (LS6 job 3600785; self-test job 3600659, 12 checks;
all 13 gates re-checked and passed). Tables in `exports/agent1/lpa_horizon/posthoc/`.
The verdict of record is unchanged and could not have been changed by these analyses.
Hazard ratios are per 144 nmol/L higher predicted Lp(a), full follow-up, unless stated.

**P1. Qualifying condition.** Mean predicted Lp(a) is 58.6 nmol/L in prior infarction,
56.8 in peripheral disease, 50.9 in prior ischaemic stroke and 48.8 in primary
prevention, so composition does differ as suspected. Adjusting for it changes little.
Four-point composite: pooled 1.19 (1.11 to 1.27), stratified 1.17 (1.09 to 1.25); stage
ratio 0.58 pooled, 0.52 (0.29 to 0.74) stratified; the verdict rule applied to the
stratified model still returns PARTIAL_PERSISTENCE. Three-point composite: 1.11 pooled,
1.10 (1.02 to 1.18) stratified, stage ratio 0.36 (0.06 to 0.66). Components, stratified:
revascularization 1.44 (1.31 to 1.58), infarction 1.26 (1.14 to 1.40), cardiovascular
death 1.07 (0.95 to 1.20), ischaemic stroke 0.86 (0.75 to 0.99). The inverse stroke
association is therefore NOT a composition artefact: pooled 0.84, stratified 0.86, and
it is seen within prior infarction (0.85, 0.70 to 1.02), prior stroke (0.93, 0.75 to
1.16) and peripheral disease (0.74, 0.49 to 1.12). It is absent at five years (0.97,
0.76 to 1.23) and opposite to primary prevention (1.14, 1.06 to 1.23). It is unexplained;
candidates are index-event selection on stroke-specific risk factors, displacement by
earlier coronary events and treatment, and chance (one of six components, nominal P about
0.03). Revascularization is raised within every condition (1.36 to 1.70). Cardiovascular
death is flat within prior infarction (1.02, 0.89 to 1.17) and raised within prior stroke
(1.30, 1.05 to 1.61).

**P2. Calibration against LDL-C.** One SD of the LDL-C score corresponds to 0.135 mmol/L
of measured LDL-C in cohort A (77% on statins) and 0.229 in cohort C; with LDL-C divided
by 0.7 in statin users, 0.236 and 0.305. Method check in cohort C: 144 nmol/L of Lp(a)
has the same association with the four-point composite as 0.78 mmol/L of LDL-C (0.63 to
0.98; 1.04 on the statin-corrected scale). The published equivalences are 101.5 mg/dL
(Burgess 2018, doi 10.1001/jamacardio.2018.1470) and 65.7 mg/dL (Lamina and Kronenberg
2019, doi 10.1001/jamacardio.2019.1041) of Lp(a) per 1 mmol/L, which for 144 nmol/L is
roughly 0.6 to 1.1 mmol/L depending on the paper and the mass-to-molar conversion: the
check passes. In cohort A the equivalence is 0.36 mmol/L (0.19 to 0.73) on measured
LDL-C and 0.63 (0.33 to 1.28) statin-corrected. Predicted hazard ratio for a trial
lowering Lp(a) by 144 nmol/L, four-point composite with any revascularization: 0.91
(0.83 to 0.95) on measured LDL-C and 0.86 (0.73 to 0.92) statin-corrected; stratified
0.91 and 0.85. Against the trial's 1.02 (0.90 to 1.15): P 0.13 and 0.04 pooled, 0.14 and
0.07 stratified. The uncalibrated lifelong estimate was 0.84 (P 0.006 against the trial).
Five-year calibrations are uninformative (the LDL-C score's five-year effect in A is too
imprecise; Fieller intervals very wide or unbounded).

**Reading.** Adjusting for the qualifying condition does not change the result. Allowing
for treatment duration moves the genetic prediction from 0.84 to roughly 0.85 to 0.91,
which the trial's interval partly covers; the trial's endpoint also excluded elective
revascularization, where the genetic effect is largest. The remaining disagreement with
the trial is borderline and depends on how statin-treated LDL-C is scaled.

## Further post hoc analyses (2026-11-20, after P1 and P2 had been read)

Declared before any of them is run. Everything above had been seen, including the
unexplained inverse stroke association and the dependence of the calibration on how
statin-treated LDL-C is scaled. None of these can change the verdict of record. Code:
`lpa_horizon_posthoc2.R`, which sources `lpa_horizon.R` unchanged and takes the function
definitions of `lpa_horizon_posthoc.R` without running it.

**P3. Is the inverse stroke association selection or displacement?** In cohort A, on
later ischaemic stroke. (a) The selection-bias budget of the pre-registered run, repeated
for each endpoint with the factor effects estimated on that endpoint and with two factors
added that the budget lacked: atrial fibrillation before baseline (first-occurrence field
131350) and baseline systolic pressure (field 4080, mean of the automated readings).
(b) The stroke model adjusted for those factors, pooled and stratified, compared with the
unadjusted model in the same participants. (c) Stroke as a first event: follow-up
censored at the first later infarction or revascularization, after which treatment
changes. Reading fixed now: if (b) or (c) moves the hazard ratio (0.86 stratified) to
within its own interval of 1 and the budget in (a) has the same sign, that mechanism is a
sufficient explanation; if neither moves it, the association stays unexplained. Limit:
the budget covers measured factors only.

**P4. Stronger statin correction.** P2 divided LDL-C by 0.7 in statin users. Repeat the
calibration with 0.6 (the 40% relative effect estimated by Chong et al. 2026, doi
10.1371/journal.pgen.1012366) and 0.5 (high-intensity treatment), and with a scale that
needs no divisor: the LDL-C difference per SD of score among participants not taking a
statin, in the primary-prevention cohort and in cohort A. Expected direction, stated
now: a stronger correction raises the LDL-C equivalence of 144 nmol/L of Lp(a), lowers
the predicted trial hazard ratio and so widens the disagreement with the trial. The
non-user scale is biased the other way, because not being treated selects for lower
LDL-C.

**P5. Is the association linear into the trial's range?** The trial enrolled at 70 mg/dL
or more (about 149 nmol/L); 7.5% of cohort A is predicted to be at or above 150 nmol/L.
(a) No outcome: by band of predicted Lp(a) (below 50, 50 to 99, 100 to 149, 150 or more
nmol/L), the share whose measured Lp(a) is at or above 150 nmol/L or above the assay's
range. (b) Hazard ratio of each band against the lowest, and the slope implied by the top
band (its log hazard ratio scaled to 144 nmol/L by the difference in mean predicted
Lp(a)). (c) A linear spline with one knot at 150 nmol/L: slope below, slope above, and
the test of their difference; knot at 100 nmol/L as a sensitivity analysis. (d) A
quadratic term. Cohorts: A pooled, A stratified and primary prevention standardized;
four-point and three-point composites; full follow-up. The slope above the knot is the
uncalibrated genetic estimate inside the trial's range and is compared with the trial as
before. Limit: about 1,300 participants of cohort A are above the knot, so that slope
will be imprecise.

**P6. Component-by-component comparison with the trial (descriptive).** Added after the
trial's Table 2 was re-read beside results already in hand, which showed the disagreement
concentrated in infarction and revascularization, so this is a description of something
already seen and not a test. For each end point the trial reports (four-component and
three-component composites, cardiovascular death, infarction, stroke, urgent coronary
revascularization) the pre-registered lifelong estimate in cohort A, per 144 nmol/L
lower, is set beside the trial's hazard ratio, with a z statistic for the difference. No
model is fitted. Limits: lifelong against 4.7 years; the trial counts any non-fatal
stroke and urgent revascularization only, this study ischaemic stroke and any
revascularization.

**External replication** is not attempted: it needs individual-level genotypes and
follow-up after disease onset in a second cohort, which is not held on LS6.

## Outcome of the further post hoc analyses

Run 2026-11-20 at commit `eae5689` (LS6 job 3600852; self-test job 3600851, 20 checks;
all 13 gates re-checked and passed). Tables in `exports/agent1/lpa_horizon/posthoc2/`.
The verdict of record is unchanged.

**P3. Stroke.** Atrial fibrillation before baseline: 9.1% of cohort A, 1.2% of cohort C.
(a) Unmodelled measured factors imply a bias of -0.016 in the stroke log hazard ratio
(diabetes -0.008, systolic pressure -0.008), against an observed -0.15: the same sign,
about a tenth of the size. (b) Adjusting for the factors: 0.83 to 0.85 (stratified, same
15,558 participants), no material change. (c) Stroke as a first event, 1,116 of the
1,267 strokes: 0.90 (0.77 to 1.04) stratified, 0.87 (0.75 to 1.02) pooled. By the reading
fixed in advance this counts as a sufficient explanation, since the interval now includes
1 and the budget has the same sign. Stated plainly, the point estimate moved only from
0.86 to 0.90; strokes that follow a coronary event account for roughly a third of the
inverse association, measured selection for a tenth, and the remainder cannot be told
from chance. In primary prevention censoring made no difference (1.14 and 1.15).

**P4. Statin correction.** LDL-C per SD of score in cohort A: 0.135 mmol/L as measured,
0.236 with divisor 0.7, 0.292 with 0.6, 0.371 with 0.5; among non-users, 0.297 in primary
prevention and 0.238 in cohort A. Divisor 0.6 reproduces the untreated primary-prevention
slope. Predicted trial hazard ratio (pooled model, full follow-up): 0.91 measured, 0.86
with 0.7, 0.82 (0.67 to 0.90) with 0.6, 0.78 (0.61 to 0.88) with 0.5; against the trial,
P 0.13, 0.04, 0.03 and 0.02 (stratified: 0.14, 0.07, 0.05, 0.04). The direction is the
one stated in advance: a stronger correction widens the disagreement with the trial.

**P5. Linearity.** (a) In the band predicted at 150 nmol/L or more (1,274 of cohort A,
mean predicted 199 nmol/L, 419 four-point events) 80% have measured Lp(a) at or above
150 nmol/L or above the assay range, so the band is what it claims to be. (b) Hazard
ratios against the lowest band, cohort A: 1.11, 1.19, 1.24; primary prevention: 1.16,
1.24, 1.47. (c) Spline with a knot at 150 nmol/L, log hazard ratio per 144 nmol/L higher:
primary prevention 0.298 below and 0.319 above (P for change 0.81), linear throughout;
cohort A 0.214 below and 0.013 above (SE 0.128; P for change 0.19). Per 144 nmol/L lower
within the trial's range, cohort A: 0.99 (0.77 to 1.27), against the trial P 0.82; knot
at 100 nmol/L 0.94 (0.79 to 1.12); three-point composite 1.02 (0.77 to 1.36). (d)
Quadratic term in cohort A P 0.13. Stratified models agree. So inside the range the trial
treated, the genetic estimate in established disease is compatible with no effect and
equally with the overall estimate of 0.84: non-linearity is suggested by the point
estimates and not established. Not declared in advance and so only noted: the trial moved
participants from about 229 to 65 nmol/L, across the knot, so the slope above the knot
alone is not the contrast the trial tested.

**P6. Components beside the trial.** Per 144 nmol/L lower, cohort A against the trial:
cardiovascular death 0.91 against 0.93 (P for difference 0.91); stroke 1.19 against 1.06
(0.43); infarction 0.78 against 0.99 (0.02); revascularization 0.67 (any) against 1.04
(urgent) (0.002); three-component 0.90 against 1.01 (0.16); four-component 0.84 against
1.02 (0.006). The disagreement is in infarction and revascularization.

**Reading.** Neither disease stage nor treatment duration, on its own, reconciles the
genetic estimate with the trial. What remains open is whether the association in
established disease weakens at the concentrations the trial treated; these data have 419
events there and cannot settle it.

Rounding corrected on 2026-11-20 in three places above (0.9146 is 0.91, not 0.92, for the
measured-scale prediction and for cardiovascular death); no estimate changed.

## Matched-endpoint analysis (2026-11-21, declared before any estimate below exists)

Post hoc. It cannot change the verdict of record. Code: `lpa_horizon_posthoc3.R`, which
sources `lpa_horizon.R` unchanged and takes the function definitions of the two post hoc
scripts and of `lpa_forecast.R` without running them.

**Why.** An outside review of the manuscript (21 November) asked for like-with-like
endpoints: the trial counted urgent coronary revascularization leading to hospitalization,
this study any coronary revascularization, which is where the genetic association is
largest. Record-level hospital episodes with the method of admission are not held.

**What had been seen when this was written.** Everything above. The forecast outputs,
including the three-point composite at the trial's concentrations: 0.89 (0.82 to 0.97) with
one slope and 0.94 (0.82 to 1.08) with two. A counts-only job (`urgent_feas_counts.R`, LS6
job 3602802) that tabulated no outcome against Lp(a) or genotype: in cohort A, 1,962 first
revascularizations, 786 participants with an acute-linked procedure and 1,180 with one that
is not; composites of 3,997 (three-point), 4,790 (four-point with any revascularization)
and 4,078 (four-point with acute-linked revascularization only). So the matched composite
is the three-point composite plus 81 events, and its estimate is known in advance to be
close to the three-point one. Hand arithmetic from printed numbers: the three-point
prediction with one slope against the trial's primary result gives P about 0.07. Not seen:
any estimate for acute-linked revascularization, for other revascularization, or for the
matched composite.

**Definitions.** Acute coronary codes: ICD-10 I20.0, I21, I22 and I24 in hospital data
(first date of each distinct code). A revascularization (OPCS-4 K40 to K46, K49, K50, K75)
is *acute-linked* if such a code was first recorded from 7 days before the procedure to one
day after it; windows of 14 and 30 days are sensitivity analyses. It is *other* if no such
code falls within 30 days. The *matched composite* is cardiovascular death, myocardial
infarction, ischaemic stroke or acute-linked revascularization (7 days). Time zero,
censoring, exposures, covariates and cohorts are those of the pre-registered analysis.

**Analyses.**

- M1, the main one. Matched composite in cohort A: hazard ratio per 144 nmol/L lower, and
  the prediction for the trial at its post-baseline geometric means (229.4 and 65.2 nmol/L)
  with one slope and with a second slope above 150 nmol/L, intervals from 20,000 draws
  (seed 42), each compared with the trial's primary result (1.02, 0.90 to 1.15) by a
  two-sided z test.
- M2. Stage contrast for the matched composite, cohort A against cohort C standardized,
  with the pre-registered decision rule applied to it.
- M3. Acute-linked and other revascularization in cohort A: a hazard ratio for each, and
  their difference from one stacked Cox model (a stratum per event type, every covariate
  free to differ by type, variance clustered on participant). The prediction for
  acute-linked revascularization at the trial's concentrations is set beside the trial's
  urgent coronary revascularization (1.04, 0.80 to 1.36).
- M4, sensitivity. Windows of 14 and 30 days; five-year horizon; cohort C standardized and
  cohort D; cohort A stratified on qualifying condition; predictions rescaled to trial
  duration through LDL-C.

**Gates.** The thirteen gates of the pre-registered run. The first revascularization date
rebuilt from the dated codes equals the endpoint of record for every participant of
cohorts A, C and D. The event counts equal those of job 3602802 (A: 786 acute-linked, 1,180
other, 4,078 matched composite). The Lp(a) coefficients of record are reproduced to 1e-9
(four-point 0.174533790354549, three-point 0.101635091628793).

**Reading fixed now.**

- M1. If the one-slope prediction for the matched composite does not differ from the
  trial's primary result at two-sided P below 0.05, the manuscript withdraws the statement
  that one slope is contradicted by the trial: on a matched endpoint the comparison is
  inconclusive, and the earlier conflict (P 0.003) is attributed to revascularizations the
  trial did not count. If it does differ at P below 0.05, the conflict is robust to
  matching the endpoint and the text says so. Either way the any-revascularization estimate
  stays the result of record for the pre-registered composite.
- M2. Reported with its verdict; the verdict of record does not change.
- M3. If the two revascularization estimates differ at two-sided P below 0.05, report that
  the association is stronger for the type with the larger estimate; otherwise report no
  evidence of a difference. No claim about urgency itself: the proxy is proximity to a
  newly recorded diagnosis, not an adjudicated admission.

**Limits fixed now.** The proxy undercounts: 32% of cohort A already carry an acute
coronary code before baseline, and a repeat event under an unchanged code is invisible, so
some urgent procedures are classed as other. The truth lies between the matched composite
and the any-revascularization composite, which are reported side by side. Diagnosis dates
are episode dates. M1 is close to the three-point result by construction, which had been
seen.

**Expectation, written before the run.** Matched composite in cohort A 0.89 to 0.91 per
144 nmol/L lower. Prediction with one slope 0.87 to 0.90, P against the trial between 0.05
and 0.15; with two slopes 0.92 to 0.96, P above 0.3. Stage ratio for the matched composite
near 0.4. Acute-linked revascularization 0.75 to 0.85, other revascularization 0.60 to
0.70, and I think it more likely than not that they differ at P below 0.05.

## Outcome of the matched-endpoint analysis

Run 2026-11-21 at commit `08b8a16` (LS6 job 3602926; self-test job 3602910, 21 checks).
The thirteen gates and the five gates of this analysis passed: counts equal to job 3602802,
the four-point coefficient reproduced exactly and the three-point one to 1e-15. Tables in
`exports/agent1/lpa_horizon/matched_endpoint/`. Hazard ratios are per 144 nmol/L lower
predicted Lp(a), full follow-up, cohort A pooled unless stated.

**M1.** Matched composite (4,078 events): 0.88 (0.82 to 0.95); stratified on qualifying
condition 0.89 (0.83 to 0.96); five-year horizon 0.84 (0.75 to 0.95). Prediction for the
trial with one slope 0.87 (0.80 to 0.94), against its primary result P 0.031; with two
slopes 0.92 (0.80 to 1.05), P 0.27. Windows of 14 and 30 days give the same numbers (P
0.031 and 0.032). Rescaled to trial duration the intervals are too wide to read (0.79,
0.40 to 0.93). For reference, the pre-registered composite gives 0.82 (P 0.003) and 0.88
(P 0.09), and the three-point composite against the trial's primary result 0.89 (P 0.07)
and 0.94 (P 0.38).

**Reading, by the rule fixed in advance.** The one-slope prediction differs from the trial
at P below 0.05, so the conflict is robust to matching the endpoint, and the manuscript
keeps it. It is weaker than on the pre-registered composite (P 0.03 against 0.003), and
with two slopes the comparison stays inconclusive.

**M2.** Stage contrast for the matched composite: cohort C standardized 0.78 (0.75 to
0.81), ratio of log hazard ratios 0.49 (0.20 to 0.79), P for difference 0.002; the rule
returns PARTIAL_PERSISTENCE, as for the composite of record.

**M3.** Acute-linked revascularization (786 events) 0.65 (0.56 to 0.76); other
revascularization (1,180 events) 0.65 (0.58 to 0.74); difference of log hazard ratios
0.001 (SE 0.095), P 0.99. No evidence of a difference. Prediction for acute-linked
revascularization at the trial's concentrations: 0.61 (0.52 to 0.73) with one slope and
0.65 (0.50 to 0.84) with two, against the trial's urgent coronary revascularization 1.04
(0.80 to 1.36), P 0.001 and 0.012. In cohort C the association is somewhat weaker for
acute-linked procedures (0.69 against 0.59, P 0.03); in cohort D 0.80 (0.63 to 1.01)
against 0.75 (0.63 to 0.91), P 0.70. Not declared and only noted: for the LDL-C score in
cohort A the two do differ (0.96 per SD lower for acute-linked, 0.88 for other, P 0.04).

**Cohort D.** Matched composite 0.96 (0.87 to 1.07), the same as its three-point
composite. The association of 0.88 (0.80 to 0.96) reported for cohort D on the
pre-registered composite is therefore carried by revascularization; on the matched
composite cohort D does not show one. Its interval includes the cohort A estimate.

**Against the written expectation.** Wrong on the two points that mattered. I expected the
one-slope prediction to lose significance on the matched composite (P 0.05 to 0.15); it is
0.03. I expected the association to be weaker for acute-linked than for other
revascularization (0.75 to 0.85 against 0.60 to 0.70); they are the same, 0.65. Right on
the rest: the matched composite itself (0.88 against an expected 0.89 to 0.91, just
outside), the two-slope prediction (0.92, P 0.27) and the stage ratio (0.49).

**What it means.** Counting any revascularization, where the trial counted urgent ones,
explains part of the disagreement and not all of it. The stand-in undercounts urgent
procedures and is not adjudicated, so this is a bound, not a like-for-like comparison:
the truth lies between the two composites, and on both the one-slope prediction is below
what the trial observed.

**Wording corrected after a second outside review (2026-11-21).** The stand-in is a proxy
for the trial's endpoint, not a match, and the manuscript calls it the trial-proxy
endpoint. The statements above, and in the limits fixed before the run, that the truth
lies between the two composites are withdrawn: the two estimates show how sensitive the
comparison is to the endpoint definition and do not establish bounds. No estimate changed.
