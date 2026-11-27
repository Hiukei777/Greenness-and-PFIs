# Residential greenness and pulmonary fungal infections

R code that reproduces the tables and figures of the manuscript examining
associations between residential greenness (NDVI and EVI at 250, 500 and
1000 m; mean and maximum) and pulmonary fungal infections (PFIs), using a
1:3 propensity-score-matched case–control design.

## Repository layout

```
.
├── run_all.R                 # run the whole pipeline
├── R/
│   ├── 00_setup.R            # libraries, input data, helper functions
│   ├── 01_table_s4_models.R  # crude / adjusted logistic models
│   ├── 02_psm_table1.R       # propensity-score matching, Table 1
│   ├── 03_figure1_balance.R  # maps, PS density, balance, correlation
│   ├── 04_table3_conditional_logistic.R
│   ├── 05_figure2_rcs.R      # restricted cubic splines, knot selection
│   ├── 06_figure3_subgroup.R # exposure–response by subgroup
│   ├── 07_table_s12_quartiles.R
│   ├── 08_sensitivity.R      # IPTW, CAPS, caliper, ratio, optimal matching
│   └── 09_evalue.R           # E-values
├── data/                     # input data (not distributed) - see data/README.md
└── Subgroup/                 # subgroup-matched datasets - see Subgroup/README.md
```

## Scripts → outputs

| Script | Analysis | Output |
| --- | --- | --- |
| `01_table_s4_models.R` | Crude and multivariable logistic regression, unmatched sample | `Table S4_crude.csv`, `Table S4_adjust.csv` |
| `02_psm_table1.R` | PS model + 1:3 nearest-neighbour matching (caliper 0.2 SD) | `Table1_tabUnmatched.csv`, `Table1_tabMatched.csv` |
| `03_figure1_balance.R` | Study-area maps, PS density, jitter and Love plots, correlation matrix | `Figure 1.pdf`, Figures S2–S4 |
| `04_table3_conditional_logistic.R` | Conditional logistic regression, overall and by subtype | `Table3_results_main.csv`, `TableS9-10.csv` |
| `05_figure2_rcs.R` | Restricted cubic splines, 3/4/5 knots | `TableS3_AIC_BIC.csv`, `Figure 2 RCS.pdf` |
| `06_figure3_subgroup.R` | Exposure–response curves by subgroup | `Figure 3 Subgroup.pdf` |
| `07_table_s12_quartiles.R` | Quartile ORs by subgroup | `TableS12_quartiles.xlsx`, `TableS12_quartiles_long.csv` |
| `08_sensitivity.R` | Unmatched, IPTW, CAPS, caliper 0.1, 1:1 and 1:5, optimal matching | `TableS4.csv`, `TableS4_iptw_bal.csv`, `TableS5_CAPS.csv`, `TableS6_0.1.csv`, `TableS7_Ratio_1.csv`, `TableS7_Ratio_5.csv`, `TableS8_Optimal.csv` |
| `09_evalue.R` | E-values for the main quartile contrasts | `TableS11_Evalues.csv` |

All outputs are written to the repository root and are ignored by git.

## Requirements

- R ≥ 4.2 (developed and run on R 4.5 / 4.6)
- Packages: `tidyverse`, `dplyr`, `readr`, `openxlsx`, `tableone`, `pastecs`,
  `MatchIt`, `WeightIt`, `cobalt`, `survey`, `survival`, `rms`, `sf`,
  `ggspatial`, `rnaturalearth`, `rnaturalearthdata`, `patchwork`, `corrplot`,
  `RColorBrewer`, `interactionR`, `EValue`

```r
install.packages(c("tidyverse", "openxlsx", "tableone", "pastecs", "MatchIt",
                   "WeightIt", "cobalt", "survey", "rms", "sf", "ggspatial",
                   "rnaturalearth", "rnaturalearthdata", "patchwork",
                   "corrplot", "RColorBrewer", "EValue"))
# interactionR is distributed on GitHub
# devtools::install_github("epi-zen/interactionR")
```

## Data

Individual-level records are **not** included in this repository. Place the
de-identified analysis dataset at `data/data.csv` and see
[`data/README.md`](data/README.md) for the expected variables.

`03_figure1_balance.R` additionally needs a provincial boundary file,
`china.geojson`, in the repository root. The boundary data are used under
standard map approval number GS(2024)0650 and are not redistributed here.

`06_figure3_subgroup.R` and `07_table_s12_quartiles.R` read the subgroup-matched
datasets in `Subgroup/`; see [`Subgroup/README.md`](Subgroup/README.md).

## Quick start

```bash
git clone <repository-url>
cd <repository>

# 1. input data
mkdir -p data
cp /path/to/data.csv data/data.csv

# 2. run everything (figures and tables appear in the repository root)
Rscript run_all.R

# ...or run a single step
Rscript R/05_figure2_rcs.R
```

All scripts must be run from the repository root: relative paths are resolved
against it, and every script sources the prerequisite steps it needs.

## Notes on the code

The analysis code is reproduced from the script used for the manuscript. Only
the following mechanical changes were made when the script was split into
steps:

1. absolute working-directory and input paths were replaced by the relative
   path `data/data.csv`; outputs are written to the repository root;
2. one output path was case-corrected (`subgroup/` → `Subgroup/`) so that it
   resolves on case-sensitive file systems;
3. the monolithic script was split into `R/00_setup.R` … `R/09_evalue.R`,
   with a short header on each file and `source()` guards for prerequisites.

No modelling code, covariate set, matching parameter, seed or plotting
parameter was altered.

### Revision of November 2026

The analysis script was revised after the first version of this repository was
prepared. Both changes are carried over verbatim:

1. `R/02_psm_table1.R` — all categorical covariates (`sex`, `occupation_4`,
   `alcohol_hx`, `smoking_hx`, `copd`, `diabetes`, `tuberculosis`, `malignant`,
   `haematological`, `immunosuppression`, `seasons`) are declared as factors
   before the propensity-score model and matching, instead of only
   `occupation_4` and `seasons`. For binary 0/1 variables this leaves the
   propensity score and the matched sets unchanged.
2. `R/06_figure3_subgroup.R` — the pre-existing-disease subgroup is now matched
   1:1 (exact on `pre_existing`) instead of 1:3. The other five subgroup
   datasets keep 1:3 matching.

The main analysis (1:3 matching, Tables 1–3, Figures 1–2) is unaffected by
either change.

## Citation

If you use this code, please cite the manuscript (see `CITATION.cff`).

## Licence

MIT — see [`LICENSE`](LICENSE).
