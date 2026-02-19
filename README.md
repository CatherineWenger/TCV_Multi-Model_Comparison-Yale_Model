# Evaluating the Impact and Cost-Effectiveness of Typhoid Conjugate Vaccine Schedules across Diverse Settings: a Multi-Model Comparison

This repo contains the code files for the medRxiv publication.   
<!-- **[medRxiv Publication:](https://doi.org/10.1371/journal.pgph.0004341)**-->

<!-- add medRxiv citation here -->

This study compares the cost-effectiveness of various typhoid conjugate vaccine (TCV) schedules under different scenarios. A compartmental transmission dynamic model was used to model the epidemiological output for each scenario and strategy. Then, the cost-effectiveness of each strategy was calculated across the different scenarios from the government perspective using the net-monetary benefit framework. This repository contains the code for the Yale model in this multi-model comparison. Repositories for the other models are provided in the supplementary documentation of the manuscript.


### Primary Analysis (10-year time horizon): 
#### Vaccine strategies
##### 1. No vaccination
##### 2. Routine TCV vaccination at 9 months with a catch-up campaign
##### 3. Routine TCV vaccination at 2 years with a catch-up campaign
##### 4. Routine TCV vaccination at 5 years with a catch-up campaign
##### 5. Routine TCV vaccination at 9 months + booster dose at 5 years with a catch-up campaign


#### Scenarios evaluated
##### 1. Incidence setting: medium (10-<100 cases per 100,000 person-years), high (100-<500 cases per 100,000 person-years), very high (>500 cases per 100,000 person-years)
##### 2. Vaccine waning scenario: slow waning (duration of protection ~ 40 years), fast waning (duration of protection ~ 8 years)
##### 3. Cost/outcome impact setting: Africa (high case-fatality ratio/high cost setting), Asia (low case-fatality ratio/low cost setting)
##### This resulted in 12 total scenario combinations in the primary analysis


#### <ins>Secondary Analyses:</ins>
##### 1. 20-year time horizon, including a two-booster vaccine strategy (Routine TCV vaccination at 9 months + booster dose at 5 & 10 years with a catch-up campaign)
##### 2. Primary vaccine strategies without a catch-up campaign
##### 3. Youngest age of vaccination at 15 months instead of 9 months


## File descriptions
#### *Transmission_dynamic_model folder* <br>
  **run_Parallel_1000_slow.R** runs... <br>  
  **run_Parallel_1000_fast.R** runs... <br>  
  **parallelization_slow.R** runs... <br>  
  **parallelization_fast.R** runs... <br>  
  **typhoid_ode_slow_waning.R** runs... <br>  
  **typhoid_ode_fast_waning.R** runs... <br>  
  **search_grid_slow.Rdata** is... <br>  
  **search_grid_fast.Rdata** is... <br>  
#### *Process_TDM_output folder* <br>
  **Process_Results.Rmd** runs... <br>  
  **post_count_doses_by_year.R** runs... <br>  
  **post_cases_per_year.R** runs... <br>  
#### *Cost_effectiveness_model folder* <br>
  **load_data.R** runs... <br>  
  **MakeDoses.R** runs... <br>  
  **MakeCases.R** runs... <br>  
  **CEA_Parameters.Rmd** runs... <br>  
  **CEA_Model_Calculations.Rmd** runs... <br> 


All model specifics and parameters are detailed in the manuscript and supplemental files published HERE. 
