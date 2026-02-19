################################################################################
# Evaluating the Impact and Cost-Effectiveness of Typhoid Conjugate Vaccine 
#   Schedules across Diverse Settings: a Multi-Model Comparison
# 
# Code to generate typhoid cases per year across strategies and simulations
#
# Created: Catherine Wenger Jan 2025
################################################################################

cases_per_year <- function(results,list_process,arch, scen,filepath){
  ## loop through the list of data by vaccine strategy to extract cases
  for(item_name in list_process){
    item = get(item_name)
    processed = list()
    ## loop through each epidemiological simulation 
    for(k in 1:length(item)){
  
      y = item[[k]]
      pop = rowSums(y[,2:217])                        # total population
      cuml_infections = y[,245:253] + y[,254:262]     # Get cumulative infections by age group
      new_infections = apply(cuml_infections, 2, diff)# Get weekly new infections
  
      ages = factor(c("0-2","2-5","5-9","10-14","15+", "All ages"),
                    levels=c("0-2","2-5","5-9","10-14","15+", "All ages"))
  
      df <- data.frame()
      ## Loop through post-vaccine time by year to condense weeks into years
      for(zy in 1:21){
        ## Yearly new infections
        new_infections_i = new_infections[((52*(zy-1))+1):((52*zy)),]
        ## Yearly population size
        pop_i = median(pop[((52*(zy-1))+1):((52*zy))])
        ## Yearly new infections by age group 
        new_infections_yr <- c(sum(rowSums(new_infections_i[,1:3])),  # 0-<2 years
                                sum(new_infections_i[,4]),            # 2-<5 years
                                sum(new_infections_i[,5]),            # 5-<10 years
                                sum(new_infections_i[,6]),            # 10-<15 years
                                sum(new_infections_i[,7:9]),          # 15+ years
                                sum(rowSums(new_infections_i[,1:9]))) # all ages
        ## add each year's data to the end of the dataframe
        df <- rbind(df, data.frame(simulation_id = k,
                                   year = zy,
                                   ages = ages,
                                   New_cases = new_infections_yr,
                                   population = rep(pop_i,6)))
  
      }
      ## add each simulation's yearly data to a list
      df$ages = factor(df$ages, levels=c("0-2","2-5","5-9","10-14","15+", "All ages"))
      processed[[k]] <- df
  
    }
    ## Name the object and save
    obj_name = paste0(item_name, "_Cases_per_year")
    file_name = paste0(gsub("_Results", "", item_name), "_Cases_per_year", ".Rdata")
    assign(obj_name, processed)
    save(list=obj_name, file=here("~/Process_TDM_output/Output",filepath,file_name))
  }
}