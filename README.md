# Firm dynamics in the German industrial sector: an analysis of the 2022 energy shock
Author: Valentina Botero, valentina.5917@gmail.com
Master Thesis in Partial Fulfilment of the Requirements for the Degree of Master of Science in Economics
Supervisor: JProf. Adam Hal Spencer, Ph.D.
Universität Bonn

This repository contains the MATLAB replication files for a two-firm-type
Hopenhayn model calibrated to German industry. Firms are classified as low or high
energy-intensive. The package reproduces the steady states, impulse responses,
robustness exercises, crisis simulation, and selected appendix diagnostics.

In required cases, AI tool was use to improve parts of the code, this was agreed previously
with my supervisor.

## Repository contents

```text
code/                         MATLAB scripts and helper functions
data/observed/                Monthly real energy-price and real-wage paths
data/calibration/             Firm-type calibration inputs
data/baselines/               Precomputed one-type and two-type steady states
results/generated/            New outputs created when scripts are run
results/selected_final_outputs/  Final figures and tables used for sharing
```

## Requirements

- MATLAB R2021a or newer is recommended.
- No non-MathWorks add-ons are required by the included baseline workflow

## Main scripts

- `ss_one_industry.m`: one-type steady states under the selected productivity
  discretization.
- `ss_two_firms_both_paths.m`: low- and high-energy firm steady states under
  Tauchen and Rouwenhorst discretizations.
- `energyshock_experiment.m`: baseline yearly IRFs plus persistence, shock-size, composition,
  and high-minus-low robustness exercises.
- `transition_paths.m`: monthly transition paths using the observed real
  energy-price and real-wage series from January 2021 onward
- `appendix_materials.m`: selected numerical diagnostics and IRFs

## Data and reproducibility

The precomputed `.mat` files allow the main simulations to run without first
repeating the steady-state calibration. Their generating scripts and the
calibration inputs are included. See `data/README.md` for the role of each
input file and verify the redistribution terms of any underlying external data
before making the repository public.

