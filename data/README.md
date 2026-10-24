# Input data

The individual-level dataset is **not** distributed with this repository.
Place the de-identified analysis file at `data/data.csv` (comma-separated,
UTF-8, 14,670 records). It is read by `R/00_setup.R` and ignored by git.

The de-identified dataset is available from the corresponding author on
reasonable request, subject to the ethics approval under which it was collected.

## Variables expected by the code

| Column | Description |
| --- | --- |
| `X` | row index |
| `patient_id` | de-identified record identifier |
| `group_01` | outcome: 1 = PFI case, 0 = control |
| `subtype_code` | 1 = pulmonary aspergillosis, 2 = pulmonary cryptococcosis, 3 = other PFIs, NA for controls |
| `age` | age, years |
| `sex` | 0 = male, 1 = female |
| `Urbanization` | urbanization category of the residential grid cell (`"3"` = city) |
| `Urban_01` | 1 = high urbanization, 0 = otherwise |
| `occupation_group` | occupation, free text |
| `occupation_4` | 0 = white-collar, 1 = blue-collar, 2 = retired / not employed, 3 = other |
| `alcohol_hx`, `smoking_hx` | 1 = yes, 0 = no |
| `copd`, `diabetes`, `tuberculosis`, `malignant`, `haematological`, `immunosuppression` | pre-existing conditions, 1 = yes, 0 = no |
| `pre_existing` | 1 = any of the six conditions above |
| `pop` | population size / density of the residential grid cell |
| `mean_temp0` | mean temperature, °C |
| `mean_humid` | mean relative humidity, proportion (0–1); multiplied by 100 for Table 1 |
| `seasons` | season category |
| `season_01` | 1 = warm season (Apr–Sep), 0 = cold season (Oct–Mar) |
| `months` | month of admission |
| `lon_wgs`, `lat_wgs` | residential coordinates, WGS84 |
| `EVI{250,500,1000}m_{mean,max}` | EVI in a 250 / 500 / 1000 m buffer, mean and maximum |
| `NDVI{250,500,1000}m_{mean,max}` | NDVI in a 250 / 500 / 1000 m buffer, mean and maximum |
| `*_IQR` | exposure divided by its interquartile range |
| `*_Q` | exposure quartile, `Q1`–`Q4` |
| `age_01` | 1 = aged ≥ 60 years |
| `ps` | propensity score for being a case |

`data_case` and `data_control` are derived in `R/00_setup.R`; the `_IQR` and
`_Q` columns are also recomputed there, so a file containing only the raw
columns listed above is sufficient.

## Map boundaries

`R/03_figure1_balance.R` reads `china.geojson` from the repository root
(provincial boundaries, standard map approval number GS(2024)0650). The file is
not redistributed with this repository.
