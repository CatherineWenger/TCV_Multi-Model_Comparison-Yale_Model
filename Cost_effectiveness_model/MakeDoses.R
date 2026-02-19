###############################################################################
# Evaluating the Impact and Cost-Effectiveness of Typhoid Conjugate Vaccine 
#  Schedules across Diverse Settings: a Multi-Model Comparison
# 
# Code to reconfigure doses per year by scenario into array for CEA
#
# Created: Catherine Wenger Jan 2025
###############################################################################

source(here("Cost_effectiveness_model","load_data.R"))
load_doses_by_year() # run this function loaded from load_data.R

years = c(1:21) # includes 1 year of pre-vaccination

## create array for doses that will fit in the parameter array made in CEA_Parameters.Rmd
nrdosesroutine =  nrdosescampaign = nrdosesbooster = 
       array(dim = c(length(geonames),
           length(burdennames),
           length(scennames),
           length(samplenames),
           length(vacstratnames),
           length(years)),
        dimnames = list(geonames,burdennames,scennames,
                        samplenames,vacstratnames,years))

## Loop through the incidence/waning scenarios
for(a in burdennames){
 if(a == "MediumSlow"){
  burden = c("Medium")
  waning = c("Slow")
 }else if(a == "HighSlow"){
  burden = c("High")
  waning = c("Slow")
 }else if(a == "VeryHighSlow"){
  burden = c("VeryHigh")
  waning = c("Slow")
 }else if(a == "MediumFast"){
  burden = c("Medium")
  waning = c("Fast")
 }else if(a == "HighFast"){
  burden = c("High")
  waning = c("Fast")
 }else if(a == "VeryHighFast"){
  burden = c("VeryHigh")
  waning = c("Fast")
 }
 ## loop through earliest age of vaccination (9/15mo)
 for(b in scennames){
  ## loop through vaccine strategies
  ##   each if statement identifies the vaccine strategy
  ##   then it loads in the data for that scenario and names it INPUT 
  for(c in vacstratnames){
   if(b == "9mo"){
    # for no vaccine strategies, this inserts zero doses administered
    if(c == "No_vax")
    {INPUT <- lapply(1:length(samplenames), function(x) {
     data.frame(
      year = years,
      ages = "All ages",
      routine_doses = rep(0, length(years)),
      Campaign_doses = rep(0, length(years)),
      Booster_doses = rep(0, length(years))
     )
    })}
    if(c == "Routine_U2_W")
     {var_name <- paste0("Scenario_",burden,"C_",waning,"_Doses_per_year")
     INPUT = get(var_name)}
    if(c == "Routine_2_W")
     {var_name <- paste0("Scenario_",burden,"G_",waning,"_Doses_per_year")
     INPUT = get(var_name)}
    if(c == "Routine_5_W")
     {var_name <- paste0("Scenario_",burden,"I_",waning,"_Doses_per_year")
     INPUT = get(var_name)}
    if(c == "Booster_5_W")
     {var_name <- paste0("Scenario_",burden,"K_",waning,"_Doses_per_year")
     INPUT = get(var_name)}
    if(c == "Booster_5_10_W")
     {var_name <- paste0("Scenario_",burden,"O_",waning,"_Doses_per_year")
     INPUT = get(var_name)}
    if(c == "Routine_U2_WO")
     {var_name <- paste0("Scenario_",burden,"B_",waning,"_Doses_per_year")
     INPUT = get(var_name)}
    if(c == "Routine_2_WO")
     {var_name <- paste0("Scenario_",burden,"F_",waning,"_Doses_per_year")
     INPUT = get(var_name)}
    if(c == "Routine_5_WO")
     {var_name <- paste0("Scenario_",burden,"H_",waning,"_Doses_per_year")
     INPUT = get(var_name)}
    if(c == "Booster_5_WO")
     {var_name <- paste0("Scenario_",burden,"J_",waning,"_Doses_per_year")
     INPUT = get(var_name)}
    if(c == "Booster_5_10_WO")
     {var_name <- paste0("Scenario_",burden,"N_",waning,"_Doses_per_year")
     INPUT = get(var_name)}
   }else if(b == "15mo"){
    ## for no vaccine strategies, this inserts zero doses administered
    if(c == "No_vax")
    {INPUT <- lapply(1:length(samplenames), function(x) {
     data.frame(
      year = years,
      ages = "All ages",
      routine_doses = rep(0, length(years)),
      Campaign_doses = rep(0, length(years)),
      Booster_doses = rep(0, length(years))
     )
    })}
    if(c == "Routine_U2_W")
     {var_name <- paste0("Scenario_",burden,"E_",waning,"_Doses_per_year")
     INPUT = get(var_name)}
    if(c == "Routine_2_W")
     {var_name <- paste0("Scenario_",burden,"G_",waning,"_Doses_per_year")
     INPUT = get(var_name)}
    if(c == "Routine_5_W")
     {var_name <- paste0("Scenario_",burden,"I_",waning,"_Doses_per_year")
     INPUT = get(var_name)}
    if(c == "Booster_5_W")
     {var_name <- paste0("Scenario_",burden,"M_",waning,"_Doses_per_year")
     INPUT = get(var_name)}
    if(c == "Booster_5_10_W")
     {var_name <- paste0("Scenario_",burden,"Q_",waning,"_Doses_per_year")
     INPUT = get(var_name)}
    if(c == "Routine_U2_WO")
     {var_name <- paste0("Scenario_",burden,"D_",waning,"_Doses_per_year")
     INPUT = get(var_name)}
    if(c == "Routine_2_WO")
     {var_name <- paste0("Scenario_",burden,"F_",waning,"_Doses_per_year")
     INPUT = get(var_name)}
    if(c == "Routine_5_WO")
     {var_name <- paste0("Scenario_",burden,"H_",waning,"_Doses_per_year")
     INPUT = get(var_name)}
    if(c == "Booster_5_WO")
     {var_name <- paste0("Scenario_",burden,"L_",waning,"_Doses_per_year")
     INPUT = get(var_name)}
    if(c == "Booster_5_10_WO")
     {var_name <- paste0("Scenario_",burden,"P_",waning,"_Doses_per_year")
     INPUT = get(var_name)}
   }
   ## loop across the number of simulations
   for(d in samplenames){
    ## First index creates the geographical region data (input data is the same for both)
    nrdosesroutine[1,a,b,d,c,] = nrdosesroutine[2,a,b,d,c,] = filter(INPUT[[d]],ages == "All ages")$routine_doses
    nrdosescampaign[1,a,b,d,c,] = nrdosescampaign[2,a,b,d,c,] = filter(INPUT[[d]],ages == "All ages")$Campaign_doses
    nrdosesbooster[1,a,b,d,c,] = nrdosesbooster[2,a,b,d,c,] = filter(INPUT[[d]],ages == "All ages")$Booster_doses
   }
  }
 }
}
