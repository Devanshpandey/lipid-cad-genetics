# Nature-round analyses: how each result will be read

Written 2026-10-30, before any of these analyses was run. Changing a rule after seeing results must be recorded here
with the reason. Scope: the Nature-format lipid-CAD manuscript (results/final/submission_nature).

## A. ApoB + Lp(a) multivariable MR
Question: does the Lp(a) effect survive conditioning on ApoB, given that Lp(a) particles carry apoB-100?
Primary model: ApoB, HDL-C, triglycerides, Lp(a) (the published four-exposure model with ApoB replacing LDL-C, same instrument
rules). Primary estimate: Lp(a) on CAD, UK Biobank and FinnGen R11. Sensitivity: two- and three-exposure models, cis-dense Lp(a)
instrument, statin-naive arm, instruments excluding the LPA window, censored Lp(a) GWAS when available.
- SUPPORTS separate axes: Lp(a) CI excludes 1 in UK Biobank and FinnGen, and the log OR is within 25% of the LDL-C-model value.
- ATTENUATED: CI still excludes 1 in both, but the log OR falls by 25-50%. The claim stays and the attenuation goes in the abstract.
- WEAKENED: log OR falls by more than 50% or the CI includes 1 in either cohort. The title and abstract claim must change.
- INVALID: positive control fails (the LDL-C analogue must reproduce OR 1.25/1.15 on CAD). No estimate is read until fixed.

## B. Selection into statin use
Question: is the statin × score asymmetry explained by who gets statins?
Tests: logistic statin use on each score; Cox interaction additionally adjusted for measured LDL-C (raw and /0.7); IPTW including
measured LDL-C with a balance table; melanoma PRS negative control; a selection-only simulation calibrated to the observed values.
- LDL-C attenuation may be called ROBUST only if the ratio of HRs CI excludes 1 under IPTW with every weighted SMD < 0.1, and the
  selection-only simulation does not reproduce an attenuation of that size.
- Otherwise it is NOT SEPARABLE FROM SELECTION: the asymmetry leaves the abstract as a finding and stays as a described association.
- SMD >= 0.1 means uninformative. It is not reported as support for either reading.
- Negative control interaction with CI excluding 1: every interaction in the section is treated as confounded.
- Lp(a) "no modification" may be claimed only as a bounded null, with its interval stated.

## C. Weak-instrument-robust inference on the LDL-C/ApoB contrast
Question: is the LDL-C minus ApoB contrast identified?
Primary: FinnGen CAD (no sample overlap with the exposures), biochemistry arm. Secondary: UK Biobank CAD, statin-naive and NMR arms.
- SUPPORTS "not identifiable": robust confidence set for the contrast is unbounded, or includes 0 in the primary analysis.
- CONTRADICTS: the set is bounded and excludes 0 in the primary analysis and agrees in sign across arms. The "not identifiable"
  claim must then be rewritten, because the data would favour one fraction.
- An empty Anderson-Rubin set is evidence of heterogeneity or misspecification, never of a precise estimate.
- Calibration must pass on simulated collinear exposures (coverage 93-97%) before any real result is read.

## D. Tudball design on the paper's scores
Question: does the paper's statin asymmetry survive Tudball's estimators and the repeat-visit LDL-C test?
Scores: the paper's LDL-C and Lp(a) scores alongside UKB standard PRS 26250, 26227, 26223 and melanoma 26252, using the existing
port cohort. Attribution: Tudball et al. 2026 (medRxiv), code MIT-licensed; adaptation disclosed.
- Coronary arm: unadjusted, LDL-adjusted and IPTW estimators. IPTW weighted SMD >= 0.1 is uninformative, as already pre-specified
  in the port's ANALYSIS_PLAN.md section 6.
- Repeat-visit arm: a null full-sample interaction with a positive LDL0-adjusted change score supports the collider mechanism.
- It is a prevalent-user design without primary-care records, not a target-trial emulation, and must be described as such.

## Amendments

Recorded 2026-11-01, before any of these analyses was run on real data, after code review of the scripts.

**B, which IPTW arm decides.** The rule above says "every weighted SMD < 0.1" without saying which form of LDL-C. Raw LDL-C in a
user is on treatment, so the drug itself separates users from non-users and truncated weights cannot close that gap. Requiring
balance on it would make ROBUST unreachable whatever the truth. The rule is therefore: the /0.7 arm must have every counted SMD
< 0.1 and a CI below 1, with raw LDL-C not counted in that arm because it is post-treatment; a raw arm that balances must also lie
below 1, and on its own it can never license ROBUST. verdict.csv also prints the both-arms-must-balance reading beside the verdict.
This reading is less strict than the literal sentence, and it was written knowing that the Tudball port's real arms did not balance
(weighted SMD 0.18-0.26); neither arm of this analysis has been seen.

**B, simulation.** "Reproduced" stays any cell in the grid, the reading that makes ROBUST harder. The result restricted to cells
matching the observed statin main effect within 10% is printed beside it. If no cell matches, the grid has not tested the observed
scenario and ROBUST is not available.

**C, the statistic.** "Robust confidence set" means the KOH subset set (Kleibergen, with the heterogeneity variance), as pre-specified in
nat_ldl_apob_robust.R. AR, K and LC sets are reported beside it. CONTRADICTS requires the statin-naive and NMR arms to be present.
The primary set is also reported without the LPA and APOE loci. That is context and does not enter the label; if it changes the set's
class or zero-inclusion the verdict file says LOCUS-SENSITIVE.

**A, sensitivity variants.** noLPAwin_* exclude chr6:157-163 Mb from selection, not the 160.5-161.7 Mb window, because rs531596603
(chr6:160,179,747) is locus-driven and lies outside the window. The rules for reading the primary variants are unchanged.

**External outcomes, A and C.** A palindromic variant ambiguous in either source (allele frequency 0.42-0.58) is dropped. The
published run_external_mr.R tests the UKB frequency only and has not been changed; the logs here count how many instruments that affects.

**B, the cohort-flow rows of the positive control.** Recorded 2026-11-04, after jobs 3545604 and
3549233 both stopped on them and no estimate of B had been read. Two rows the script enforced as
"published" were not published values. The start of 402,119 traces to an unrelated Lp(a) diagnostic
(logs/claude_jobs/lpa_diag_3542525.out, "merged rows 402200 in-field 402119"), and the 376 "dropped"
is the residual 402,119 - 21,847 - 379,896 rather than a reported count; the manuscript states no
cohort flow at all. Both rows are now printed and exported with no target and no pass or fail. The
real start is 401,743. Everything else stays a hard stop: the four interaction and main-effect
estimates against the manuscript (to 2 dp, plus the 5-dp full-precision pair), and the counts n,
events, statin users, events in statin users and the prevalent exclusion, which come from a prior run
of the paper's own script and so verify reproduction of that run, not agreement with a printed table.

**C, the NMR arm reads the published instrument file.** Recorded 2026-11-05, AFTER the real C result
(job 3549234, verdict UNRESOLVED) had been read, at the author's instruction. In that run the NMR arm
reselected instruments from the archived per-chromosome GWAS with this script's distance-only clump and
failed its positive control (LDL-P 2.94 and LDL-C(NMR) 0.536 against the published 2.48 and 0.69), so
the arm counted as invalid. Supplementary Table oS3b was produced from
exports/agent1/nmr_gencov/nmr_mvmr_snps.csv (62 harmonised instruments, CARDIoGRAM outcome, 2026-09-16),
so the arm now reads that file. The reason is faithfulness to the published analysis, not direction.
That file carries no alleles, so FinnGen cannot be harmonised to it and the NMR arm has the CARDIoGRAM
outcome only. The control target is the export's full-precision MVMR-IVW estimate (LDL-P 2.47578,
LDL-C(NMR) 0.686485). The verdict rule is unchanged.
Written before the rerun: on the published estimates the NMR contrast LDL-C(NMR) minus LDL-P is about
-0.371 - 0.907 = -1.28, the opposite sign to the biochemistry arm's positive set. The rule counts a
secondary arm whose qhet point estimate has the opposite sign as disagreement, so this change is expected
to leave the label at UNRESOLVED, now for disagreement across arms rather than for a missing arm, and
not to produce CONTRADICTS.

**Correction to the B cohort-flow amendment of 2026-11-04.** Recorded 2026-11-07, after B's real result had been read; it changes
no verdict. That amendment said the manuscript states no cohort flow. That was wrong: a text search for "379,896" missed the LaTeX
form 379{,}896. The previous main text printed "Excluding 21,847 participants with coronary disease on or before baseline and 376
without follow-up left 379,896", and the Supplementary Table caption (oS34) printed "from 402,119 with a baseline date, a medication
report and score data". So the two rows were published values, and they do not reproduce. The rerun reaches the same 379,896
participants, 32,053 events and estimates, but starts from 401,743 with a baseline date, a medication report and score data and
drops nobody for lack of follow-up or covariates. The 376 are most likely participants without a medication report
(402,119 - 401,743), counted under the wrong label. This is not verified here. The rows stay reported without a target, and the
manuscript flow is corrected in revision r1 instead of being treated as absent.
