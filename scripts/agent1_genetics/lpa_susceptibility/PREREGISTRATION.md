# Does arterial susceptibility amplify the Lp(a) axis as it does the ApoB axis?

Written 2026-11-21, before any estimate below is computed. Code:
`lpa_susceptibility.R`, which uses the metabolome-resilience libraries
(`scripts/agent6_resilience/R/lib_*.R`) and the Lp(a) exposure model of
`scripts/agent1_genetics/lpa_horizon/lpa_horizon.R`, all unchanged.

## Question

The metabolome-resilience work found that measured ApoB and a non-lipid coronary
polygenic score more than multiply each other's hazard (product term 0.058, SE 0.015, in
181,162 statin-free participants; 0.039, P 0.011, replicated with clinical ApoB), and
randomized trials show larger relative benefit of LDL lowering at higher genetic risk
(Mega et al. 2015, doi 10.1016/S0140-6736(14)61730-X; Natarajan et al. 2017, doi
10.1161/CIRCULATIONAHA.116.024436). Is the effect of Lp(a), the second lipid axis,
amplified in the same way?

## What had been seen when this was written

All results of the metabolome-resilience pipeline (S1 to S3 and the clinical
replication), including that the genetic version of the ApoB test, lipid-region score by
non-lipid score, was null at the primary window (0.012, SE 0.019). All results of the
Lp(a)HORIZON emulation and forecasts, in which genetically predicted Lp(a) and the
non-lipid score entered the same models as main effects only. Not seen: any product of
Lp(a), measured or predicted, with any polygenic score; any product of the 167-variant
LDL-C score with the non-lipid score.

## Data and model, taken unchanged from the earlier analysis

- Cohort: all eligible participants of its stage 11 (European ancestry, genetic quality
  control, no prevalent coronary disease, follow-up beyond a two-year landmark):
  379,654 with 15,978 incident events. Outcome: myocardial infarction or its
  complications (I21 to I23, I24.1, I25.2), coronary revascularization (OPCS-4 K40 to
  K46, K49, K50, K75) or coronary death, to 31 October 2022, clock restarted at the
  landmark.
- Scores: the lipid-region (L) and non-lipid (N) parts of the UK Biobank-free coronary
  polygenic score, regions within 500 kb of ApoB lead variants, each centred and divided
  by the SD of the total score.
- Model: Cox regression with M, L, N, their squares and the products M:L, M:N and L:N,
  adjusted for age, sex, ten principal components and genotyping batch (the "full
  quadratic" of `quad_amplification`).

## Exposure

M is genotype-predicted Lp(a) from the censored-normal model on 45 *LPA* variants
(the emulation's exposure), standardized in the eligible cohort. It is uncorrelated with
the non-lipid score (0.001 in the emulation cohort).

## Primary test and decision rule

The coefficient of M:N for predicted Lp(a), two-sided at 0.05. To compare axes the
product is expressed as r = b(M:N) / b(M), the fractional change in the exposure's own
log hazard ratio per unit of N, which does not depend on the units of M. The ApoB
reference r is re-estimated from the earlier discovery fit (expected about 0.29), and
"half" means half of that value.

| Verdict | Condition | Lead-grade |
|---|---|---|
| AMPLIFIED | 95% CI of b(M:N) above 0 and r at least half of ApoB's | yes |
| AMPLIFIED_LESS_THAN_HALF | 95% CI of b(M:N) above 0 and r below half of ApoB's | no |
| REVERSED | 95% CI of b(M:N) below 0 | yes |
| NOT_AMPLIFIED_AXES_DIFFER | upper 95% limit of r below half of ApoB's, and r differs from ApoB's at P < 0.01 | yes |
| INCONCLUSIVE | anything else | no |

The axis contrast is r(Lp(a)) minus r(ApoB), with delta-method standard errors and the
two estimates treated as independent (the exposures are nearly uncorrelated; the samples
overlap). A lead-grade verdict is called robust only if the statin-free and spline
versions give a product of the same sign whose interval lies on the same side of zero.

## Secondary analyses (reported in full, no verdict)

1. The M:L product and the earlier contrast D.
2. Products beyond 4-df spline main effects.
3. Statin-free participants; each sex; score windows of 250 and 1,000 kb.
4. Measured Lp(a), with values below and above the assay range set to its limits.
5. The same model with the 167-variant LDL-C score as M: the genetic comparator for the
   ApoB axis.
6. Additive scale: relative excess risk due to interaction for the top decile of
   predicted Lp(a) with the top quintile of N.
7. One model in the earlier discovery sample holding measured ApoB and predicted Lp(a)
   together, with the difference of their r from a single covariance matrix.

## Gates (nothing is read if any fails)

The thirteen gates of the emulation (exposure model and cohort); the eligible cohort has
379,654 participants and 15,978 events, and its statin-free part 324,626 and 11,740; the
earlier lipid-region by non-lipid product is reproduced (0.01244, tolerance 1e-6); the
discovery sample has 181,162 participants and 6,540 events and the earlier ApoB product
is reproduced (0.05766, tolerance 1e-6); at least 99.9% of the eligible cohort has
predicted Lp(a); the correlation of predicted Lp(a) with N is below 0.02 in absolute
value.

## Limits fixed in advance

- The ApoB reference rests on measured ApoB, which is not randomized; its genetic
  version was null. The comparison of a measured exposure with a genetic one is
  therefore unequal, and secondary analysis 5 is there to show the genetic comparator.
- A product on the log-hazard scale can be absent while risks still add; analysis 6
  reports the additive scale.
- Predicted Lp(a) is strongly skewed, so the product may be carried by its upper tail;
  the squares in the model and the spline version address curvature, not this.
- One cohort, European ancestry, healthy-volunteer selection.
- The non-lipid against lipid distinction of the earlier work did not replicate
  independently.

## Expectation, written before the run

The product for Lp(a) will be close to zero, r between -0.10 and +0.10. Most likely
verdicts: INCONCLUSIVE or NOT_AMPLIFIED_AXES_DIFFER. Probability of AMPLIFIED about one
in five. The LDL-C score comparator will also be close to zero.

## Deviations

None at the time of writing. Any change after this file is committed is listed here.

## Outcome

Run 2026-11-21 at commit `62b1861` (LS6 job 3602201; self-test job 3602197, 16 checks).
All gates passed: the thirteen emulation gates, the recorded cohort and event counts, and
exact reproduction of the two earlier estimates (0.0124399151 and 0.0576577396). Tables
in `exports/agent1/lpa_susceptibility/`.

**Verdict: INCONCLUSIVE, not lead-grade.** Product of predicted Lp(a) with the non-lipid
score 0.0072 (SE 0.0087; 95% CI -0.0098 to 0.0241; P 0.41). As a fraction of the Lp(a)
effect, r = 0.065 (-0.094 to 0.225). The ApoB reference re-fitted in the discovery sample
gave r = 0.290 (SE 0.084), so half of it is 0.145: the Lp(a) interval reaches above that,
and the contrast between the axes is -0.225, P 0.054.

**Secondary.** Spline version 0.0081 (SE 0.0082); statin-free 0.0061 (0.0100); windows
of 250 and 1,000 kb 0.0066 and 0.0056; sex coded 0 0.0236 (0.0150) and sex coded 1
-0.0018 (0.0106). Measured Lp(a) 0.0042 (0.0090), r 0.024 (-0.079 to 0.127). The
167-variant LDL-C score, the genetic comparator for the other axis, -0.0099 (0.0099): no
amplification either. One model with both exposures in the discovery sample: r 0.286 for
ApoB and 0.037 for Lp(a), difference 0.249 (SE 0.131), P 0.058. Additive scale, top
decile of predicted Lp(a) with top quintile of the score: hazard ratios 1.12, 2.16 and
2.50, relative excess risk 0.22 (SE 0.105), multiplicative departure 0.031 (SE 0.050):
the risks multiply, so they more than add.

**A flaw in secondary analysis 1, found on reading the result.** The product of
predicted Lp(a) with the lipid-region score was -0.046 (SE 0.016; P 0.004), and D was
0.099 (SE 0.034). The lipid regions are those around ApoB lead variants, and three of
those leads lie at *LPA*, among them rs10455872. The lipid-region score therefore
contains the *LPA* locus, the same locus that defines the exposure, so this product and
D are not a test of Lp(a) against other lipid loci and are not interpreted. I had noted
this risk while scoping and did not carry it into the design. A lipid-region score with
the *LPA* locus removed would be needed, and would be post hoc.

**Against the written expectation.** r was 0.065, inside the expected -0.10 to +0.10.
The verdict is one of the two expected. The LDL-C score comparator was close to zero as
expected.

**Reading.** There is no evidence that arterial susceptibility amplifies the Lp(a) axis:
Lp(a) and non-lipid polygenic risk multiply. The data do not show that the axes differ
(P about 0.05 by two routes), and the genetic comparator for the ApoB axis shows no
amplification, so the earlier amplification remains a property of measured ApoB.

## Post hoc analysis (2026-11-21, after the outcome above had been read)

Declared before it is run. It cannot change the verdict, which rests on the product with
the non-lipid score. Code: `lpa_susceptibility_posthoc.R`, which takes the function
definitions of `lpa_susceptibility.R` without running it.

**Seen:** the product of predicted Lp(a) with the lipid-region score was -0.046 (P
0.004), and the lipid-region score contains the *LPA* locus (lipid region chr6:159,991,741
to 162,021,039 at the 500 kb window, holding three ApoB lead variants, among them
rs10455872).

**Analysis.** The part of the lipid-region score carried by variants inside that one
region is scored separately with the same PLINK2 options, weights and genotype files as
the original components, and subtracted. Three models in the eligible cohort, each the
full quadratic with the same covariates:

1. the original, to confirm the pre-registered product is reproduced;
2. the lipid-region score with the *LPA* region removed: the product of predicted Lp(a)
   with it is the corrected version of secondary analysis 1, with the corrected D;
3. the *LPA*-region component entered as its own term beside the remainder, to show
   where the original product came from.

Model 2 is repeated in statin-free participants and with measured Lp(a). Reported
without an outcome: the correlations of predicted Lp(a) with each component and the
share of the lipid-region score's variance that the *LPA* region carries.

**Reading fixed now.** If the product in model 2 has a 95% interval that includes zero
and the product with the *LPA*-region component in model 3 is negative, the original
result was the same locus on both sides of a product and nothing more. If the product in
model 2 stays negative with an interval excluding zero, Lp(a) combines less than
multiplicatively with lipid-locus polygenic risk, which would then need replication
before it is believed.

**Gates.** The eligible cohort and its events as before; the original product reproduced
to 1e-6; the *LPA*-region score present for at least 99.9% of the cohort; the two parts
sum to the original lipid-region score for every participant; exactly one lipid region
overlaps the *LPA* gene and it contains rs10455872.

## Outcome of the post hoc analysis

Run 2026-11-21 at commit `8cc1ada` (LS6 job 3602413; self-test job 3602296, 8 checks).
Gates passed: cohort and events as before, both pre-registered products reproduced to
ten decimal places, the two parts summing to the original score to 5e-16. The *LPA*
region held 1,066 score variants. Tables in
`exports/agent1/lpa_susceptibility/posthoc_lpa_removed/`.

**Components, no outcome.** The *LPA* region carries 20% of the variance of the
lipid-region score. Predicted Lp(a) correlates 0.30 with the original lipid-region score,
0.67 with its *LPA*-region part and -0.003 with the remainder.

**Products with the lipid-region score** (per SD of Lp(a) by unit of score):

| Analysis | Original | *LPA* region removed | *LPA* region as its own term |
|---|---|---|---|
| Predicted Lp(a), all eligible | -0.046 (SE 0.016), P 0.004 | -0.015 (0.017), P 0.38 | -0.119 (0.043), P 0.005 |
| Predicted Lp(a), statin-free | -0.055 (0.019), P 0.004 | -0.020 (0.020), P 0.29 | -0.110 (0.052), P 0.035 |
| Measured Lp(a) | -0.036 (0.016), P 0.030 | -0.016 (0.018), P 0.36 | 0.064 (0.044), P 0.14 |

The corrected contrast D is 0.041 (SE 0.032), against 0.099 (0.034) before. The product
with the non-lipid score is unchanged (0.0085, SE 0.0082), so the verdict stands.

**Reading, by the rule fixed in advance.** For predicted Lp(a), the analysis that was
pre-registered, the corrected product includes zero and the product with the
*LPA*-region component is negative: the original result was the same locus on both sides
of the product. There is no evidence that Lp(a) combines other than multiplicatively
with polygenic risk at the other lipid loci. One wrinkle: with measured Lp(a) the
corrected product is equally null, but the same-locus term is not negative, so the
"same locus" account is supported for predicted Lp(a) and not shown for measured.

**Note for the earlier work.** One fifth of the variance of its lipid-region score is
the *LPA* locus. Its ApoB analyses are not invalidated by this, since *LPA* is an ApoB
locus, but "lipid-region" there means LDL-pathway loci plus Lp(a).
