# Summary-level outputs

Copies of the tables the jobs wrote, with their provenance files (commit, job, date, md5 of the code). The md5 of each R script
named in a provenance file equals the md5 of the copy in this repository. No table holds individual-level data, an identifier,
or a participant or event count from 1 to 4.

| Folder | Contents |
|---|---|
| `lpa_horizon/feasibility`, `forecast_feasibility`, `urgent_feasibility` | Counts made before the analyses, with no outcome tabulated against genotype |
| `lpa_horizon/emulation` | The comparison: estimates, stage contrast, verdict, gates |
| `lpa_horizon/posthoc`, `posthoc2` | Post hoc checks: composition, stroke, LDL-C scaling, dose shape |
| `lpa_horizon/matched_endpoint` | The trial-proxy endpoint (the folder keeps the name used when it was run) |
| `lpa_horizon/forecast` | Forecasts over the grid of baselines and lowering, dose models, the check against Lp(a)HORIZON |
| `lpa_susceptibility` | The amplification test and the rerun with the *LPA* region removed |
| `nature_round/coloc_susie5` | Multi-signal colocalization: the run reported, and the first run with its log (`run1_prespecified`) |
