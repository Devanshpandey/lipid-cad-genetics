# Parts of a separate pipeline that this paper uses

| Files | Use in this paper |
|---|---|
| `prs/`, `slurm/02*`, `R/10a_prs_partition.R`, `slurm/10_synergy.slurm` | A coronary polygenic score built without UK Biobank (PRS-CS on a meta-analysis of the Million Veteran Program, CARDIoGRAMplusC4D and FinnGen) and split into the part within 500 kb of ApoB-associated loci and the rest. The comparison with Lp(a)HORIZON adjusts for the second part |
| `R/01_build_cohort.R`, `config.yaml` | The cohort the amplification test starts from. `config.yaml` is the settings file of the whole pipeline, so it names data and stages not used here |
| `R/lib_*.R`, `config_synergy.yaml` | Libraries and settings the amplification test loads; their md5 are checked before it runs |

`slurm/10_synergy.slurm` goes on to two later stages (`10_synergy.R`, `11_genetic_synergy.R`) that belong to a metabolomics analysis
not reported in this paper and not included here. What this paper needs (`prs_components.tsv`) is written before them. The measured-
ApoB comparator of the amplification test also reads two intermediate files of that analysis (`03/nmr_technical_resid.tsv` and
`04/mrs_oof.tsv`).
