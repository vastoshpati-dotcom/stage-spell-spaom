# SS-SPAOM Monte Carlo Simulation

This folder contains the final Monte Carlo simulation program and summary
outputs for the Stage–Spell Peace Agreement Opportunity Model (SS-SPAOM).

## Main program

- `02_SS_SPAOM_Revised_Simulation_FINAL.R`

## Simulation design

The main study contains 96 scenarios formed by:

- baseline design: common and unequal;
- number of peace processes: 100, 184, 300, and 500;
- mean process length: 8, 12, and 20;
- stage-spell parameter: 0, 0.25, 0.52, and 0.80.

Each scenario used 300 Monte Carlo replications.

## Main output files

- `01_scenario_grid.csv`: complete scenario design;
- `05_parameter_recovery_summary.csv`: recovery of the spell parameter;
- `06_eta_recovery_summary.csv`: recovery of origin-stage intensities;
- `07_model_performance_summary.csv`: predictive performance summaries;
- `09_paired_advantage_summary.csv`: paired predictive comparisons;
- `10_best_model_regions.csv`: regions in which each model performed best;
- `SIMULATION_INTERPRETATION.txt`: concise interpretation of the results;
- `R_session_info.txt`: software and package information.

Positive paired predictive advantages favour the SS-SPAOM.
