# Subgroup-matched datasets

`R/06_figure3_subgroup.R` and `R/07_table_s12_quartiles.R` read matched
datasets from this directory. They are **not** distributed with the repository,
for the same data-sharing reasons as `data/data.csv`.

| File | Subgroup variable | Used for |
| --- | --- | --- |
| `matchdatam_age.csv` | `age_01` (< 60 / ≥ 60 years) | Figure 3 a–b, Table S12 |
| `matchdatam_sex.csv` | `sex` (male / female) | Figure 3 c–d, Table S12 |
| `matchdatam_occu4.csv` | `occupation_4` (4 categories) | Figure 3 e–f, Table S12 |
| `matchdatam_urban01.csv` | `Urban_01` (low / high urbanization) | Figure 3 g–h, Table S12 |
| `matchdatam_season01.csv` | `season_01` (cold / warm season) | Figure 3 i–j, Table S12 |
| `matchdatam_pre_existing.csv` | `pre_existing` (yes / no) | Figure 3 k–l, Table S12 |

Each file has the same columns as `data.csv` plus the matching output of
`MatchIt` (`distance`, `weights`, `subclass`). They are produced by the
subgroup-specific propensity-score matching step: the same propensity-score
model, caliper and 1:3 ratio as the main analysis, with exact matching on the
subgroup variable. `R/06_figure3_subgroup.R` regenerates
`matchdatam_pre_existing.csv` itself (see the first block of that script), so
that file is optional; the other five are required inputs.

Provide the files with these exact names if you obtained them together with the
de-identified dataset, or regenerate them with the same matching specification
described in the Methods.
