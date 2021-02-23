# /inputs

This folder contains inputs to the VSEIR model built from scripts in ``./code/R/``

The input files are:
* ``age_group_contact_matrix.csv`` : Mean daily contacts from [Prem et al. (2020)](https://www.medrxiv.org/content/10.1101/2020.07.22.20159772v2)
* ``national_baseline.csv`` : Group level characteristics (e.g., population size, efficacy, hesitancy)
* ``pums_data_with_nasem_phases.csv`` : ACS observations labeled to a unique NASEM phase
* ``seir_model_nat_lvl_infection_starting_points_by_group_on_dec14.csv`` : Group level share susceptible, exposed, infectious, and resolved as of December 14, 2020
* ``vaccination_flow_strategy_nasem_unreserved_with_[XX]_uptake_param.csv`` : Vaccine roll-out under NASEM guidelines under parameter age-stratified hesitancy assumption