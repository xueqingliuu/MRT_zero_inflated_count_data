# MRT_zero_inflated_count_data
This repository contains the code for conducting the simulation study and data analysis for the paper `Incorporating nonparametric methods for estimating causal excursion effects in mobile health with zero-inflated count outcomes` by Xueqing Liu, Tianchen Qian, Lauren Bell, and Bibhas Chakraborty.

## Simulation Study
As described in Section 5, we have four simulation scenarios. Under each scenario, there are mainly four files:
-`estimator.R` contains the implementation of estimators EMEE, EMEE-NonP, DR-EMEE-NonP, GEE, ECE, ECE-NonP, G-estimator of Yu et al. 2023
-`data_generation.R` generates the simulated MRT data 
-`other_functions.R` contains some functions used in the above files
-`run_simulation_marginal` is the main file to perform the simulation for marginal causal effects
-`run_simulation_moderation` is the main file to perform the simulation for fully conditional causal effects

## Drink Less Data Analysis
-`estimator.R` contains the implementation of estimators EMEE, EMEE-NonP, DR-EMEE-NonP for two-arm MRTs
-`estimator_three_cate.R` contains the implementation of estimators EMEE, EMEE-NonP, DR-EMEE-NonP for three-arm MRTs
-`drinkless.R` is the main file to perform the analysis of Drink Less MRT data

