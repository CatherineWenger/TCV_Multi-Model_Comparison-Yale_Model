###############################################################################
# Evaluating the Impact and Cost-Effectiveness of Typhoid Conjugate Vaccine 
#  Schedules across Diverse Settings: a Multi-Model Comparison
# 
# Code to reconfigure cases per year by scenario into array for CEA
#
# Created: Catherine Wenger Jan 2025
###############################################################################

source(here("Cost_effectiveness_model","load_data.R"))
load_cases_by_year() # run this function loaded from load_data.R

years = c(1:21)

## create array for cases that will fit in the parameter array made in CEA_Parameters.Rmd
cases = cases0_2 = cases2_5 = cases5_10 = cases10_15 = cases15over = population = 
  array(dim = c(length(geonames),
        length(incidencenames),
        length(scennames),
        length(samplenames),
        length(vacstratnames),
        length(years)),
        dimnames = list(geonames,incidencenames,scennames,
                        samplenames,vacstratnames,years))

## Loop through the incidence/waning scenarios
##    psymp - the fraction of primary cases that are symptomatic (fit parameter)
##    cases are multiplied by psymp to get only clinical infections for the CEA
for(a in incidencenames){
 if(a == "MediumSlow"){
  psymp = 0.028
  incidence = c("Medium")
  waning = c("Slow")
 }else if(a == "HighSlow"){
  psymp = 0.112
  incidence = c("High")
  waning = c("Slow")
 }else if(a == "VeryHighSlow"){
  psymp = 0.605
  incidence = c("VeryHigh")
  waning = c("Slow")
 }else if(a == "MediumFast"){
  psymp = 0.028
  incidence = c("Medium")
  waning = c("Fast")
 }else if(a == "HighFast"){
  psymp = 0.112
  incidence = c("High")
  waning = c("Fast")
 }else if(a == "VeryHighFast"){
  psymp = 0.605
  incidence = c("VeryHigh")
  waning = c("Fast")
 }
 ## loop through earliest age of vaccination (9/15mo)
 for(b in scennames){
  ## loop through vaccine strategies
  ##   each if statement identifies the vaccine strategy
  ##   then it loads in the data for that scenario and names it INPUT 
  for(c in vacstratnames){
   if(b == "9mo"){
    if(c == "No_vax")
     {var_name <- paste0("Scenario_",incidence,"A_Slow_Results_Cases_per_year")
      INPUT = get(var_name)}
    if(c == "Routine_U2_W")
     {var_name <- paste0("Scenario_",incidence,"C_",waning,"_Results_Cases_per_year")
     INPUT = get(var_name)}
    if(c == "Routine_2_W")
     {var_name <- paste0("Scenario_",incidence,"G_",waning,"_Results_Cases_per_year")
     INPUT = get(var_name)}
    if(c == "Routine_5_W")
     {var_name <- paste0("Scenario_",incidence,"I_",waning,"_Results_Cases_per_year")
     INPUT = get(var_name)}
    if(c == "Booster_5_W")
     {var_name <- paste0("Scenario_",incidence,"K_",waning,"_Results_Cases_per_year")
     INPUT = get(var_name)}
    if(c == "Booster_5_10_W")
     {var_name <- paste0("Scenario_",incidence,"O_",waning,"_Results_Cases_per_year")
     INPUT = get(var_name)}
    if(c == "Routine_U2_WO")
     {var_name <- paste0("Scenario_",incidence,"B_",waning,"_Results_Cases_per_year")
     INPUT = get(var_name)}
    if(c == "Routine_2_WO")
     {var_name <- paste0("Scenario_",incidence,"F_",waning,"_Results_Cases_per_year")
     INPUT = get(var_name)}
    if(c == "Routine_5_WO")
     {var_name <- paste0("Scenario_",incidence,"H_",waning,"_Results_Cases_per_year")
     INPUT = get(var_name)}
    if(c == "Booster_5_WO")
     {var_name <- paste0("Scenario_",incidence,"J_",waning,"_Results_Cases_per_year")
     INPUT = get(var_name)}
    if(c == "Booster_5_10_WO")
     {var_name <- paste0("Scenario_",incidence,"N_",waning,"_Results_Cases_per_year")
     INPUT = get(var_name)}
   }else if(b == "15mo"){
    if(c == "No_vax")
     {var_name <- paste0("Scenario_",incidence,"A_Slow_Results_Cases_per_year")
     INPUT = get(var_name)}
    if(c == "Routine_U2_W")
     {var_name <- paste0("Scenario_",incidence,"E_",waning,"_Results_Cases_per_year")
     INPUT = get(var_name)}
    if(c == "Routine_2_W")
     {var_name <- paste0("Scenario_",incidence,"G_",waning,"_Results_Cases_per_year")
     INPUT = get(var_name)}
    if(c == "Routine_5_W")
     {var_name <- paste0("Scenario_",incidence,"I_",waning,"_Results_Cases_per_year")
     INPUT = get(var_name)}
    if(c == "Booster_5_W")
     {var_name <- paste0("Scenario_",incidence,"M_",waning,"_Results_Cases_per_year")
     INPUT = get(var_name)}
    if(c == "Booster_5_10_W")
     {var_name <- paste0("Scenario_",incidence,"Q_",waning,"_Results_Cases_per_year")
     INPUT = get(var_name)}
    if(c == "Routine_U2_WO")
     {var_name <- paste0("Scenario_",incidence,"D_",waning,"_Results_Cases_per_year")
     INPUT = get(var_name)}
    if(c == "Routine_2_WO")
     {var_name <- paste0("Scenario_",incidence,"F_",waning,"_Results_Cases_per_year")
     INPUT = get(var_name)}
    if(c == "Routine_5_WO")
     {var_name <- paste0("Scenario_",incidence,"H_",waning,"_Results_Cases_per_year")
     INPUT = get(var_name)}
    if(c == "Booster_5_WO")
     {var_name <- paste0("Scenario_",incidence,"L_",waning,"_Results_Cases_per_year")
     INPUT = get(var_name)}
    if(c == "Booster_5_10_WO")
     {var_name <- paste0("Scenario_",incidence,"P_",waning,"_Results_Cases_per_year")
     INPUT = get(var_name)}
   }
   ## loop across the number of simulations
   for(d in samplenames){
    ## ages: 0-2  2-5   5-9   10-14  15+   All ages
    ## First index creates the geographical region data (input data is the same for both)
    cases[1,a,b,d,c,] = dplyr::filter(INPUT[[d]],ages == "All ages")$New_cases*psymp
    cases[2,a,b,d,c,] = dplyr::filter(INPUT[[d]],ages == "All ages")$New_cases*psymp
    cases0_2[1,a,b,d,c,] = dplyr::filter(INPUT[[d]],ages == "0-2")$New_cases*psymp
    cases0_2[2,a,b,d,c,] = dplyr::filter(INPUT[[d]],ages == "0-2")$New_cases*psymp
    cases2_5[1,a,b,d,c,] = dplyr::filter(INPUT[[d]],ages == "2-5")$New_cases*psymp
    cases2_5[2,a,b,d,c,] = dplyr::filter(INPUT[[d]],ages == "2-5")$New_cases*psymp
    cases5_10[1,a,b,d,c,] = dplyr::filter(INPUT[[d]],ages == "5-9")$New_cases*psymp
    cases5_10[2,a,b,d,c,] = dplyr::filter(INPUT[[d]],ages == "5-9")$New_cases*psymp
    cases10_15[1,a,b,d,c,] = dplyr::filter(INPUT[[d]],ages == "10-14")$New_cases*psymp
    cases10_15[2,a,b,d,c,] = dplyr::filter(INPUT[[d]],ages == "10-14")$New_cases*psymp
    cases15over[1,a,b,d,c,] = dplyr::filter(INPUT[[d]],ages == "15+")$New_cases*psymp
    cases15over[2,a,b,d,c,] = dplyr::filter(INPUT[[d]],ages == "15+")$New_cases*psymp
    population[1,a,b,d,c,] = dplyr::filter(INPUT[[d]],ages == "All ages")$population
    population[2,a,b,d,c,] = dplyr::filter(INPUT[[d]],ages == "All ages")$population
   }
  }
 }
}


