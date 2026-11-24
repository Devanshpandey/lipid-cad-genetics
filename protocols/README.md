# Protocol documents and their history

Four documents record what was fixed before an analysis was run, and what was added afterwards. They stay beside the code they
govern and are reproduced verbatim.

| Document | Governs |
|---|---|
| `scripts/agent1_genetics/lpa_horizon/PREREGISTERED_READOUT.md` | The genetic comparison with Lp(a)HORIZON: protocol, gates, decision rule, outcome, then every post hoc analysis with the date it was declared and what it found, including two wrong predictions |
| `scripts/agent1_genetics/lpa_horizon/FORECAST_PREREGISTRATION.md` | The forecasts for ACCLAIM-Lp(a) and OCEAN(a)-Outcomes and the rules for judging them |
| `scripts/agent1_genetics/lpa_susceptibility/PREREGISTRATION.md` | The test of polygenic amplification of the Lp(a) axis |
| `scripts/agent1_genetics/nature_round/PREREGISTERED_READOUT.md` | Analyses added in review (Supplementary Note 3). The ApoB analysis had been read before this document was committed |

## What the history shows and does not show

`history/` holds, for each commit that touched one of these documents, the commit hash, the dates, the message and the change to
the document. `history/commits.tsv` lists the scripts committed at the same time with their md5. `history/after_last_commit.diff`
holds any text added after the last commit; it is empty in v2.0.0. The outcome sections of the three documents, written on 20
and 21 November 2026, were committed on 24 November, shortly before this release (`history/12_*`).

These commits were made in a private working repository. The hashes quoted in the paper (238b4c3, c2f7eca, 08b8a16) are from
there and do not appear in this repository's own log. **A commit is not a public timestamp:** its date is set by the author's
computer and nobody outside could see it. The forecasts count as made in advance only from the date on which this repository's
release, a preprint or a registration first made them public.

## Words used in these documents

They were written during the work and are not edited to match the paper.

| In the documents | In the paper |
|---|---|
| pre-registered | pre-specified, or specified before outcome modelling (nothing was publicly registered) |
| emulation | genetic comparison with the trial |
| matched endpoint, matched composite | trial-proxy endpoint. The word "matched" overstated it; a dated note in the readout says so and withdraws a sentence that read the two composites as bounds |
| cohorts A, C, D | established disease at baseline; before a first event; from 30 days after a first event in follow-up |

Job numbers refer to the Lonestar6 cluster.
