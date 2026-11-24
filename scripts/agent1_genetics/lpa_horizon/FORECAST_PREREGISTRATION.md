# Forecasts for Lp(a) outcome trials that have not reported

Written 2026-11-20, before any forecast below is computed. Companion to
`PREREGISTERED_READOUT.md` (the genetic emulation of Lp(a)HORIZON). Code:
`lpa_forecast.R`, which sources `lpa_horizon.R` unchanged and takes function definitions
from the two post hoc scripts without running them.

A local commit is not a public timestamp. For these forecasts to count as made in
advance they must appear in a preprint or a public registration before the trials
report; that step is the authors' and is not done by this file.

## Question

Two trials of near-complete Lp(a) lowering are running. What do the genetic estimates in
UK Biobank predict for them, stated now, and what result would contradict the
prediction?

## What had been seen when this was written

Everything in `PREREGISTERED_READOUT.md`: the stage contrast (Lp(a) association after
disease onset 58% of that before), the components, the dose spline in the
established-disease cohort (slope above 150 nmol/L 0.99 per 144 nmol/L lower, 0.77 to
1.27), and the LDL-C scaling. A counts-only table for the cohorts defined below (job
3600899), with no outcome tabulated against Lp(a). Rough hand arithmetic from printed
estimates: about 0.73 (linear) and 0.80 (spline) for an OCEAN(a)-like contrast on the
four-point composite, and about 0.88 (0.77 to 0.99) for the Lp(a)HORIZON contrast under
the spline. Not seen: any estimate in the cohorts that include prior angioplasty or
bypass grafting, in the at-risk cohort, on the coronary endpoints, or any spline outside
the two emulation cohorts.

## The trials (ClinicalTrials.gov, read 2026-11-20)

| | OCEAN(a)-Outcomes, NCT05581303 | ACCLAIM-Lp(a), NCT06292013 |
|---|---|---|
| Drug | olpasiran | lepodisiran |
| Participants | 7,297 | 17,300 |
| Lp(a) at entry | 200 nmol/L or more | 175 nmol/L or more |
| Population | prior infarction and/or coronary angioplasty with a further risk factor | (a) established disease with an event or revascularization; (b) age 55 or more, at risk of a first event |
| Primary end point | coronary death, infarction or urgent coronary revascularization | cardiovascular death, infarction, ischaemic stroke or urgent coronary revascularization |
| Pre-specified by stratum | no | yes: established disease, at risk of a first event, Lp(a) 200 nmol/L or more |
| Primary completion | March 2028 | March 2029 |

Lowering in phase 2: olpasiran about 97% from a median 260 nmol/L (O'Donoghue et al.
2022, doi 10.1056/NEJMoa2211023); lepodisiran about 94% from a median 254 nmol/L (Nissen
et al. 2025, doi 10.1056/NEJMoa2415818). The baseline distributions of the phase 3 trials
have not been published as far as a PubMed search on 2026-11-20 found.

## UK Biobank analogues

Exposure, covariates, censoring and event definitions are those of the emulation. Two
end points are added: a coronary composite (death with underlying cause I20 to I25,
infarction, any coronary revascularization) and its hard version without
revascularization. Prior angioplasty is OPCS-4 K49, K50 or K75 on or before baseline or
self-reported operation 1070; prior bypass grafting is K40 to K46 or operation 1095.

| Trial stratum | Cohort | n | Events | End points |
|---|---|---|---|---|
| Lp(a)HORIZON (check) | prior infarction, ischaemic stroke or peripheral disease | 16,955 | 4,790 four-point | four-point, three-point |
| OCEAN(a)-Outcomes | prior infarction and/or angioplasty | 12,128 | 3,455 coronary composite | coronary composite (primary), hard coronary, three-point |
| ACCLAIM established | prior event, angioplasty or bypass grafting | 20,140 | 6,006 four-point | four-point (primary), three-point |
| ACCLAIM at risk, primary | no vascular disease, age 55 or more, with diabetes or at least two of diabetes, hypertension, current smoking | counted at run time | | four-point (primary), three-point |
| ACCLAIM at risk, sensitivity 1 | no vascular disease, age 55 or more | 228,783 | 20,071 four-point | four-point, three-point |
| ACCLAIM at risk, sensitivity 2 | the primary definition plus participants aged 55 or more with vascular disease recorded but no event or revascularization | counted at run time | | four-point |
| ACCLAIM at risk, sensitivity 3 | all primary prevention (the emulation's cohort C) | 371,735 | 24,875 four-point | four-point |

The trial's "combination of high-risk factors" is not spelled out in the registry, so
the at-risk analogue is a judgement; the primary definition is fixed here and the others
are reported beside it. Cohorts are not weighted to one another.

## Models

One Cox model per cohort and end point, with predicted Lp(a), the LDL-C score, the
non-lipid coronary score, age, sex and ten principal components, in two forms:

- **linear**: log hazard proportional to predicted Lp(a);
- **spline**: a second slope above 150 nmol/L.

For a contrast from level L0 (untreated) to L1 (treated), the log hazard ratio is the
linear slope times (L0 - L1), or under the spline the two slopes applied to the parts of
the fall above and below 150 nmol/L. Each is reported on two scales:

- **lifelong**: the genetic estimate as it stands;
- **scaled to trial duration**: through the LDL-C equivalent in the same cohort (LDL-C in
  statin users divided by 0.6) and the rate ratio of 0.78 per mmol/L of the statin
  trials.

Intervals are the 2.5th and 97.5th percentiles of 20,000 draws from the fitted
coefficient distribution (seed 42); the uncertainty of the 0.78 and of the LDL-C slope
is ignored.

## Scenarios for the Lp(a) levels

- Lp(a)HORIZON: the published time-averaged post-baseline levels, 229.4 and 65.2 nmol/L.
- OCEAN(a)-Outcomes and ACCLAIM-Lp(a): baseline 225, 250, 275 or 300 nmol/L and
  time-averaged lowering of 80%, 90% or 97%. Central scenario, an assumption: 275 nmol/L
  and 90% for OCEAN(a)-Outcomes, 250 nmol/L and 90% for ACCLAIM-Lp(a), the same in both
  of its strata. The forecast to compare with each trial is the grid cell nearest to what
  it reports.

## Forecasts and how each will be judged

1. **Primary: the stage contrast inside ACCLAIM-Lp(a).** R, the log hazard ratio in the
   established-disease stratum divided by that in the at-risk stratum, on the four-point
   composite, linear model, lifelong scale. Under that model R does not depend on the
   Lp(a) levels if they are the same in the two strata. The forecast is the value of R
   and its interval; it is contradicted if the trial's at-risk stratum does not show the
   larger relative benefit, or if the observed ratio lies outside the interval after
   allowing for the trial's own uncertainty. If neither stratum shows benefit, R is
   undefined and the absolute forecasts are contradicted instead.
2. **Secondary: the absolute effect in OCEAN(a)-Outcomes.** Hazard ratio for the coronary
   composite in the prior-infarction-or-angioplasty cohort, both models and both scales,
   over the grid. Contradicted if the trial's 95% interval excludes the forecast interval
   of every model and scale at the nearest grid cell.
3. **Secondary: absolute effects in the two ACCLAIM-Lp(a) strata**, likewise.
4. **Check: Lp(a)HORIZON in retrospect.** The same machinery at the published levels,
   against the observed 1.02 (0.90 to 1.15) and, for the three-point composite, 1.01
   (0.88 to 1.16). This is not a forecast; it shows how the method fared where the answer
   is known.

## Gates (no forecast is read if any fails)

The thirteen gates of the emulation, and: the OCEAN(a)-like cohort has 12,128
participants and 3,455 coronary-composite events; the ACCLAIM established cohort 20,140
and 6,006 four-point events; the age-55-or-more cohort without vascular disease 228,783
and 20,071.

## Limits fixed in advance

- **Extrapolation.** Only 309 participants of the OCEAN(a)-like cohort are predicted at
  200 nmol/L or more, and 1,241 of the ACCLAIM established cohort at 175 nmol/L or more.
  Every forecast applies slopes estimated mostly below the trials' range.
- The linear model already over-predicted Lp(a)HORIZON.
- The trials count urgent revascularization; this study any revascularization, where the
  genetic association is largest. The three-point and hard coronary versions are given
  for that reason.
- A lifelong genetic difference is not a drug; the scaling assumes Lp(a) behaves like
  LDL-C in that respect.
- Established-disease cohorts are survivors of their first event.
- The phase 3 Lp(a) levels are assumptions until the trials publish them.

## Expectation, written before the run

R between 0.5 and 0.7. OCEAN(a)-Outcomes coronary composite about 0.6 to 0.7 under the
linear model on the lifelong scale and nearer 0.75 to 0.8 under the spline. The
Lp(a)HORIZON check under the spline about 0.85 to 0.9, not significantly different from
the observed result.

## Deviations

None at the time of writing. Any change after this file is committed is listed here.

## Outcome

Run 2026-11-20 at commit `c2f7eca` (LS6 job 3600935; self-test job 3600930, 14 checks;
the thirteen emulation gates and the six forecast gates passed). Tables in
`exports/agent1/lpa_horizon/forecast/`. Hazard ratios are for each trial's lowering at
its central scenario, lifelong scale, with 95% intervals.

**Check, Lp(a)HORIZON (229.4 to 65.2 nmol/L).** Four-point composite: one slope 0.82
(0.76 to 0.88), against the observed 1.02 P 0.003; two slopes 0.88 (0.77 to 0.99), P
0.09. Scaled to trial duration: 0.80 (P 0.02) and 0.86 (P 0.13). Three-point: 0.89 and
0.94, neither contradicted (P 0.13 and 0.47). The one-slope model fails where the answer
is known; the two-slope model is not rejected.

**Primary forecast, ACCLAIM-Lp(a) stage contrast.** R = 0.69 (0.40 to 1.16) with the
primary at-risk definition (20,065 participants, 3,261 events); probability that the
at-risk stratum benefits more, 92%. The interval includes 1: the forecast gives a
direction, not a clear separation. Broader at-risk definitions: 0.59 (0.38 to 0.82) for
all aged 55 or more without vascular disease, 0.57 (0.37 to 0.79) for all primary
prevention, 0.67 (0.41 to 1.02) adding recorded vascular disease without an event.
Three-point composite 0.44 (0.14 to 0.95).

**Absolute forecasts.** ACCLAIM-Lp(a) at risk: 0.68 (0.59 to 0.78) one slope, 0.67 (0.53
to 0.86) two slopes. ACCLAIM-Lp(a) established: 0.77 (0.70 to 0.84) and 0.84 (0.74 to
0.97). OCEAN(a)-Outcomes coronary composite: 0.71 (0.62 to 0.80) and 0.81 (0.65 to 1.02);
coronary death or infarction 0.75 (0.65 to 0.87) and 0.88 (0.67 to 1.15); three-point
0.86 (0.75 to 0.98) and 1.06 (0.83 to 1.36).

**Slopes above 150 nmol/L in the wider cohorts** (not a declared test; printed with the
models): change of slope -0.24 (SE 0.13, P 0.07) in the ACCLAIM established cohort and
-0.28 (SE 0.17, P 0.10) for the coronary composite in the OCEAN(a)-like cohort; none in
any at-risk cohort (P 0.4 to 0.99). These cohorts contain the emulation's
established-disease cohort and are not independent of it.

**Against the written expectation.** R 0.69, inside 0.5 to 0.7. OCEAN(a)-Outcomes 0.71
with one slope, at the edge of the expected 0.6 to 0.7, and 0.81 with two, just above the
expected 0.75 to 0.8. The Lp(a)HORIZON check under two slopes 0.88 and not significantly
different, as expected.

**What a null result would now mean.** Stated while scoping, and wrong: that a null
OCEAN(a)-Outcomes would contradict both models. The two-slope interval reaches 1.02 for
its coronary composite and 1.15 without revascularization, so a null result there would
contradict only the one-slope model, which Lp(a)HORIZON has already contradicted. The
forecast that a null result would clearly contradict is the at-risk stratum of
ACCLAIM-Lp(a): 0.67 to 0.68, upper limit 0.86 under either model.

**Reliability of the scaled forecasts.** Where the LDL-C score's association with the
end point is weak (three-point composites in established disease) the scaled intervals
are very wide and should not be used.
