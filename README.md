# Firm Dynamics in Germany Following an Energy-Price Shock

This repository contains the MATLAB replication files for a two-firm-type
Hopenhayn model of German industry. Firms are classified as low or high
energy-intensive. The package reproduces the steady states, impulse responses,
robustness exercises, crisis simulation, and selected appendix diagnostics.

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
- No non-MathWorks add-ons are required by the included baseline workflow.
- A complete run of the IRF and robustness script can take several minutes.

## Quick start

1. Open MATLAB and set the Current Folder to the repository root.
2. Run `setup_project.m` once per MATLAB session.
3. Use the precomputed steady states in `data/baselines/`, or recalibrate them
   using the first two scripts below.
4. Run the desired analysis script.

```matlab
setup_project

% Optional: rebuild one-type and two-type steady states.
run(fullfile('code', 'ss_one_industry.m'))
run(fullfile('code', 'ss_two_firms_both_paths.m'))

% Baseline IRFs and robustness exercises.
run(fullfile('code', 'test.m'))

% Monthly crisis simulation using observed input-price paths.
run(fullfile('code', 'Crisis_Simulation.m'))

% Appendix tables and figures. Run test.m first to create the required IRFs.
run(fullfile('code', 'appendix_materials.m'))
```

New files are written to `results/generated/`. The files under
`results/selected_final_outputs/` have stable descriptive names and are not
overwritten by the scripts.

## Main scripts

- `ss_one_industry.m`: one-type steady states under the selected productivity
  discretization.
- `ss_two_firms_both_paths.m`: low- and high-energy firm steady states under
  Tauchen and Rouwenhorst discretizations.
- `test.m`: baseline yearly IRFs plus persistence, shock-size, composition,
  and high-minus-low robustness exercises.
- `Crisis_Simulation.m`: monthly transition paths using the observed real
  energy-price and real-wage series from January 2021 onward.
- `appendix_materials.m`: selected numerical diagnostics and full omitted IRFs.

## Data and reproducibility

The precomputed `.mat` files allow the main simulations to run without first
repeating the steady-state calibration. Their generating scripts and the
calibration inputs are included. See `data/README.md` for the role of each
input file and verify the redistribution terms of any underlying external data
before making the repository public.

## Citation

Complete the author surname and GitHub username in `CITATION.cff` before
publishing. GitHub will then provide a **Cite this repository** option.

Suggested thesis text:

```latex
Replication files are available at
\url{https://github.com/YOUR-USERNAME/two-type-energy-shock-matlab}.
```

## License

No reuse license has been selected yet. The included `LICENSE` file reserves
all rights until the thesis author deliberately replaces it with the preferred
license, such as MIT for code intended for reuse.
