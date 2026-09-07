# Data

The raw dataset used in the paper is **not shared in this repository** due to institutional data-sharing constraints. It is available from the corresponding author on reasonable request.

To run the pipeline on your own comparably structured data, place an `.xlsx` file in this folder (or point `data_path` in `R/pipeline_final_12june.R` at it) with the following columns.

## Expected schema

One row per animal per observation.

| Column       | Type    | Description |
|--------------|---------|-------------|
| `animal_no`  | integer | Unique animal identifier. Used for animal-aware train/test splitting. |
| `Status`     | factor  | One of `healthy`, `CM`, `SCM`. Ground-truth label. |
| `DOB`        | date    | Date of birth. |
| `DOR`        | date    | Date of last calving (used to compute DIM). |
| `DOC`        | date    | Date of observation / imaging. |
| `T1`–`T4`    | numeric | Mean surface temperature (°C) of udder quarters 1–4. |
| `U1`–`U4`    | numeric | Mean surface temperature (°C) of teat quarters 1–4. |
| `T1max`–`T4max` | numeric | Maximum temperature within each udder quarter ROI. |
| `U1max`–`U4max` | numeric | Maximum temperature within each teat ROI. |
| `AT`         | numeric | Ambient temperature (°C) at the time of imaging. |
| `RH`         | numeric | Ambient relative humidity (%). |
| `Parity`     | integer | Lactation number. |

Ground-truth labels in the paper were derived from somatic cell count (SCC) and California Mastitis Test (CMT), not from bacteriological culture. See the manuscript Methods for the cut-offs used.

## Features engineered from these columns

The `engineer_features()` function in the pipeline derives, per row:

- Inter-quarter differentials for T and U (each quarter minus the mean of the others).
- Inter-quarter ranges (max − min) for T, U, Tmax, Umax.
- Inter-quarter standard deviations for T, U, Tmax, Umax.
- Maximum single-quarter differential.
- Lactation stage from `DOC − DOR`: `Early` (≤ 90 d), `Mid` (91–180 d), `Late` (> 180 d).
- Binary and three-class labels for the four modelling tasks.

Rows with any missing required column are dropped before modelling.
