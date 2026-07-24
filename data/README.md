# Data files

## `observed/`

- `Pet.csv`: monthly real energy-price index used in the crisis simulation.
- `wt.csv`: monthly real-wage index used in the crisis simulation.

Both files contain a date column and one numeric series. The crisis script
reads the complete files and selects the analysis period internally.

## `calibration/`

- `gamma_classification_rules.csv`: energy shares and firm-count shares under
  alternative sector classification rules.
- `gamma_aggregate_BC_C.csv`: aggregate energy-share calibration summary.
- `sector_classification_median_split.csv`: sector classification underlying
  the baseline low/high split.

## `baselines/`

The four `.mat` files are precomputed one-type and two-type steady states under
Tauchen and Rouwenhorst discretizations. They are included so the main IRF and
crisis scripts can run immediately. They can be regenerated with
`ss_one_industry.m` and `ss_two_firms_both_paths.m`.

Before publishing, add the original data-source citations and confirm that the
observed and calibration inputs may be redistributed publicly.
