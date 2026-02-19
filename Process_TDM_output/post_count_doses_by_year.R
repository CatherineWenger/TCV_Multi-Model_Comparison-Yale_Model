################################################################################
# Evaluating the Impact and Cost-Effectiveness of Typhoid Conjugate Vaccine 
#   Schedules across Diverse Settings: a Multi-Model Comparison
# 
# Code to generate doses administered per year across strategies and simulations
#
# Created: Catherine Wenger Jan 2025
################################################################################

# Doses administered per year


count_doses_by_year <- function(results,list_process,arch, scen,filepath){
  ## remove no vaccination strategy since doses will be zero
  results <- results[-1]
  list_process <- list_process[-1]
  list_process1 <- c()
  ## loop through the list of data by vaccine strategy to extract cases
  for (item_name in list_process) {
  
    item = get(item_name)
    # Initialize a matrix to hold doses data for each age group
    list_doses = list()
  
    ## loop through each epidemiological simulation 
    for (l in 1:length(item)) {
  
      y = item[[l]]
      cdoses = apply(y[,218:226], 2, diff) # Campaign doses
      bdoses = apply(y[,227:235], 2, diff) # Booster doses
      rdoses = apply(y[,236:244], 2, diff) # Routine doses
  
      df <- data.frame()
      ## Loop through post-vaccine time by year to condense weeks into years
      for(z in 1:21){
        ## convert weekly doses into yearly totals
        cdoses_i = cdoses[((52*(z-1))+1):((52*z)),] # Campaign doses
        bdoses_i = bdoses[((52*(z-1))+1):((52*z)),] # Booster doses
        rdoses_i = rdoses[((52*(z-1))+1):((52*z)),] # Routine doses
  
        ## Yearly campaign doses administered by age group 
        camp_doses <- c(sum(rowSums(cdoses_i[,1:3])),      # 0-<2 years
                             sum(cdoses_i[,4]),            # 2-<5 years
                             sum(cdoses_i[,5]),            # 5-<10 years
                             sum(cdoses_i[,6]),            # 10-<15 years
                             sum(cdoses_i[,7:9]),          # 15+ years
                             sum(rowSums(cdoses_i[,1:9]))) # all ages
        ## Yearly booster doses administered by age group 
        booster_doses <- c(sum(rowSums(bdoses_i[,1:3])), # 0-<2 years
                        sum(bdoses_i[,4]),               # 2-<5 years
                        sum(bdoses_i[,5]),               # 5-<10 years
                        sum(bdoses_i[,6]),               # 10-<15 years
                        sum(bdoses_i[,7:9]),             # 15+ years
                        sum(rowSums(bdoses_i[,1:9])))    # all ages
        ## Yearly routine doses administered by age group 
        routine_doses <- c(sum(rowSums(rdoses_i[,1:3])), # 0-<2 years
                        sum(rdoses_i[,4]),               # 2-<5 years
                        sum(rdoses_i[,5]),               # 5-<10 years
                        sum(rdoses_i[,6]),               # 10-<15 years
                        sum(rdoses_i[,7:9]),             # 15+ years
                        sum(rowSums(rdoses_i[,1:9])))    # all ages
        ## add each year's data to the end of the dataframe
        df <- rbind(df, data.frame(simulation_id = l,
                                    year = z,
                                    ages = c("0-1","2-4","5-9","10-14","15+","All ages"),
                                    Campaign_doses = camp_doses,
                                    Booster_doses = booster_doses,
                                    routine_doses = routine_doses))
  
      }
  
      ## add each simulation's yearly data to a list
      df$ages = factor(df$ages, levels=c("0-1","2-4","5-9","10-14","15+","All ages"))
      list_doses[[l]] <- df
      ## Name the object and save
      obj_name = paste0(gsub("_Results", "", item_name), "_Doses_per_year")
      file_name = paste0(gsub("_Results", "", item_name), "_Doses_per_year", ".Rdata")
      assign(obj_name, list_doses)
      save(list=obj_name, file=here("~/Process_TDM_output/Output",filepath,file_name))
    }
  }
}
