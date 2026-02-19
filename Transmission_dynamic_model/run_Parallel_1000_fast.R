###############################################################################
# Evaluating the Impact and Cost-Effectiveness of Typhoid Conjugate Vaccine Schedules across Diverse Settings: a Multi-Model Comparison 
# 
# Execution code to run all vaccine and incidence simulations under fast-waning assumptions
#
# Created: Catherine Wenger Jan 2025
###############################################################################

library(here)
library(tidyverse)
library(doSNOW)
library(foreach)

##Vaccination Scenes:::
#   A - No vaccination
#   B - Routine @ 9mo without catch-up
#   C - Routine @ 9mo with catch-up
#   D - Routine @ 15mo without catch-up
#   E - Routine @ 15mo with catch-up
#   F - Routine @ 2yrs without catch-up
#   G - Routine @ 2yrs with catch-up
#   H - Routine @ 5yrs without catch-up
#   I - Routine @ 5yrs with catch-up
#   J - Routine @ 9mo without catch-up + 5yrs Booster
#   K - Routine @ 9mo with catch-up + 5yrs Booster
#   L - Routine @ 15mo without catch-up + 5yrs Booster
#   M - Routine @ 15mo with catch-up + 5yrs Booster
#   N - Routine @ 9mo without catch-up + 5 & 10yrs Booster
#   O - Routine @ 9mo with catch-up + 5 & 10yrs Booster
#   P - Routine @ 15mo without catch-up + 5 & 10yrs Booster
#   Q - Routine @ 15mo with catch-up + 5 & 10yrs Booster

t0=5200;              # Length of burn-in period (weeks)
tmod = (52*20) + 1    # Model for 20 years post vaccine introduction
tspan=t0+tmod;        # Total modeled time horizon

number_runs = 1000    # number of simulations

## create output save directory
dir.create("~/Transmission_dynamic_model/Output", recursive = TRUE, showWarnings = FALSE)

# set up parallelization for running simulations
num_cores <- parallel::detectCores() - 1
cl <- makeCluster(num_cores, type = "SOCK")
on.exit(stopCluster(cl))
print("Registering parallel backend")
registerDoSNOW(cl)
print("Parallel backend registered")
clusterExport(cl, c( "tspan"))

for(i in 1:3){
  # Run for the Medium Incidence setting (i=1)
  if(i == 1){   
  load(paste0("~/Transmission_dynamic_model/","search_grid_fast.Rdata"))                                          # load 1000 samples for stochastic parameters
    scenes = c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q')           # run all vaccine scenarios

    source(here("parallelization_fast.R"))                                                         # load in "run_simulations_parallel" function
    log_file_path <- c(paste0("~/Transmission_dynamic_model/Output/","scenario_run_times_log_full_run.txt"))  # Create log file to track progress of model simulations
    
    # start log message and export parameters to the workers/cores
    num_runs = number_runs

    log_message <- paste0("------------New Runs----------------: ", Sys.Date())
    write(log_message, file = log_file_path, append = TRUE)
    start_time_overall <- Sys.time()
    clusterExport(cl, c("scenes", "num_runs", "run_simulations_parallel", "log_file_path", "here"))
    
    # Parallel loop through each scenario for 100 
    all_results <- foreach(scene = scenes,.export = c("tspan","scenes", "here"), .packages = c("here", "dplyr")) %dopar% {  
      tryCatch({

      al = 9 # number of age groups
      arch = c("Medium") # incidence archetype for these simulations

      # Initialize a list and dataframe to store the simulation results
      scenario_results <- list()
      params_df <- data.frame()

      for (p in 1:num_runs) {
        # Record the start time
        start_time <- Sys.time()

        # Run the simulation function
        sim_output <- run_simulations_parallel(p, scene, arch)

        # Record the end time and calculate the duration
        end_time <- Sys.time()
        duration <- end_time - start_time

        # Store the results
        scenario_results[[p]] <- sim_output$results

        # Combine the stochastic (non-fixed) parameters into one dataframe
        params_df <- rbind(params_df, c(Simulation_ID = p, sim_output$non_fixed_parms))

        # Print the scenario, iteration, and duration
        log_message <- paste0("Scenario ", arch, scene, " Fast | ", p, " of ", num_runs,  " | Duration: ", round(duration, 2), " seconds", " | Time Stamp: ", Sys.time())

        # Print the log message
        print(log_message)

        # Append the log message to the log file
        write(log_message, file = log_file_path, append = TRUE)
      }
      
      # Assign the simulation results and parameters to a list
      return(list(results = scenario_results, params = params_df))
      # Explicitly call garbage collection to free up memory
      gc()
  })
    }
    # Finish writing the log file
    end_time_overall <- Sys.time()
    duration_overall <- end_time_overall - start_time_overall
    log_message <- paste0("------------Done---------------- | Total Time: ", duration_overall)
    write(log_message, file = log_file_path, append = TRUE)
    
    # Save all vaccine scenario simulations together
    for(idx in seq_along(all_results)) {
      results <- all_results[[idx]]$results
      params <- all_results[[idx]]$params
      
      scene <- scenes[idx]
      results_var_name <- paste0("Scenario_Medium", scene, "_Fast_Results")
      params_var_name <- paste0("Scenario_Medium", scene, "_Fast_Params")
      
      # Construct the file name dynamically for the results and parameters
      file_name_results <- paste0("~/Transmission_dynamic_model/Output/Scenario_Medium", scene, "_Fast_Results.Rdata")
      file_name_params <- paste0("~/Transmission_dynamic_model/Output/Scenario_Medium", scene, "_Fast_Params.Rdata")
      
      # Save the results and parameters separately
      save(list = "results", file = file_name_results)
      save(list = "params", file = file_name_params)
    }
  # Run for the High Incidence setting (i=2)  
  } else if (i == 2){
      load(paste0("~/Transmission_dynamic_model/","search_grid_fast.Rdata"))                                        # load 1000 samples for stochastic parameters
      scenes = c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q')           # run all vaccine scenarios
      
      source(here("parallelization_fast.R"))                                                         # load in "run_simulations_parallel" function
      log_file_path <- c(paste0("~/Transmission_dynamic_model/Output/","scenario_run_times_log_full_run.txt"))  # Create log file to track progress of model simulations
      
      num_runs = number_runs
      
      # start log message and export parameters to the workers/cores
      log_message <- paste0("------------New Runs----------------: ", Sys.Date())
      write(log_message, file = log_file_path, append = TRUE)
      start_time_overall <- Sys.time()
      clusterExport(cl, c("scenes", "num_runs", "run_simulations_parallel", "log_file_path", "here"))
      
      # Parallel loop through each scenario for 100 
      all_results <- foreach(scene = scenes,.export = c("tspan","scenes", "here"), .packages = c("here", "dplyr")) %dopar% {  
        tryCatch({
          
          al = 9 # number of age groups
          arch = c("High") # incidence archetype for these simulations
          
          # Initialize a list and dataframe to store the simulation results
          scenario_results <- list()
          params_df <- data.frame()
          
          for (p in 1:num_runs) {
            # Record the start time
            start_time <- Sys.time()
            
            # Run the simulation function
            sim_output <- run_simulations_parallel(p, scene, arch)
            
            # Record the end time and calculate the duration
            end_time <- Sys.time()
            duration <- end_time - start_time
            
            # Store the results
            scenario_results[[p]] <- sim_output$results
            
            # Combine the stochastic (non-fixed) parameters into one dataframe
            params_df <- rbind(params_df, c(Simulation_ID = p, sim_output$non_fixed_parms))
            
            # Print the scenario, iteration, and duration
            log_message <- paste0("Scenario ", arch, scene, " Fast | ", p, " of ", num_runs,  " | Duration: ", round(duration, 2), " seconds", " | Time Stamp: ", Sys.time())
            
            # Print the log message
            print(log_message)
            
            # Append the log message to the log file
            write(log_message, file = log_file_path, append = TRUE)
          }
          
          # Assign the simulation results and parameters to a list
          return(list(results = scenario_results, params = params_df))
          # Explicitly call garbage collection to free up memory
          gc()
        })
      }
      # Finish writing the log file
      end_time_overall <- Sys.time()
      duration_overall <- end_time_overall - start_time_overall
      log_message <- paste0("------------Done---------------- | Total Time: ", duration_overall)
      write(log_message, file = log_file_path, append = TRUE)
      
      # Save all vaccine scenario simulations together
      for(idx in seq_along(all_results)) {
        results <- all_results[[idx]]$results
        params <- all_results[[idx]]$params
        
        scene <- scenes[idx]
        results_var_name <- paste0("Scenario_Medium", scene, "_Fast_Results")
        params_var_name <- paste0("Scenario_Medium", scene, "_Fast_Params")
        
        # Construct the file name dynamically for the results and parameters
        file_name_results <- paste0("~/Transmission_dynamic_model/Output/Scenario_Medium", scene, "_Fast_Results.Rdata")
        file_name_params <- paste0("~/Transmission_dynamic_model/Output/Scenario_Medium", scene, "_Fast_Params.Rdata")
        
        # Save the results and parameters separately
        save(list = "results", file = file_name_results)
        save(list = "params", file = file_name_params)
      }
      # Run for the Very High Incidence setting (i=2)  
  } else if (i == 3){
    load(paste0("~/Transmission_dynamic_model/","search_grid_fast.Rdata"))                                          # load 1000 samples for stochastic parameters
    scenes = c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q')           # run all vaccine scenarios
    
    source(here("parallelization_fast.R"))                                                         # load in "run_simulations_parallel" function
    log_file_path <- c(paste0("~/Transmission_dynamic_model/Output/","scenario_run_times_log_full_run.txt"))  # Create log file to track progress of model simulations
    
    # start log message and export parameters to the workers/cores
    num_runs = number_runs
    
    log_message <- paste0("------------New Runs----------------: ", Sys.Date())
    write(log_message, file = log_file_path, append = TRUE)
    start_time_overall <- Sys.time()
    clusterExport(cl, c("scenes", "num_runs", "run_simulations_parallel", "log_file_path", "here"))
    
    # Parallel loop through each scenario for 100 
    all_results <- foreach(scene = scenes,.export = c("tspan","scenes", "here"), .packages = c("here", "dplyr")) %dopar% {  
      tryCatch({
        
        al = 9 # number of age groups
        arch = c("VeryHigh") # incidence archetype for these simulations
        
        # Initialize a list and dataframe to store the simulation results
        scenario_results <- list()
        params_df <- data.frame()
        
        for (p in 1:num_runs) {
          # Record the start time
          start_time <- Sys.time()
          
          # Run the simulation function
          sim_output <- run_simulations_parallel(p, scene, arch)
          
          # Record the end time and calculate the duration
          end_time <- Sys.time()
          duration <- end_time - start_time
          
          # Store the results
          scenario_results[[p]] <- sim_output$results
          
          # Combine the stochastic (non-fixed) parameters into one dataframe
          params_df <- rbind(params_df, c(Simulation_ID = p, sim_output$non_fixed_parms))
          
          # Print the scenario, iteration, and duration
          log_message <- paste0("Scenario ", arch, scene, " Fast | ", p, " of ", num_runs,  " | Duration: ", round(duration, 2), " seconds", " | Time Stamp: ", Sys.time())
          
          # Print the log message
          print(log_message)
          
          # Append the log message to the log file
          write(log_message, file = log_file_path, append = TRUE)
        }
        
        # Assign the simulation results and parameters to a list
        return(list(results = scenario_results, params = params_df))
        # Explicitly call garbage collection to free up memory
        gc()
      })
    }
    # Finish writing the log file
    end_time_overall <- Sys.time()
    duration_overall <- end_time_overall - start_time_overall
    log_message <- paste0("------------Done---------------- | Total Time: ", duration_overall)
    write(log_message, file = log_file_path, append = TRUE)
    
    # Save all vaccine scenario simulations together
    for(idx in seq_along(all_results)) {
      results <- all_results[[idx]]$results
      params <- all_results[[idx]]$params
      
      scene <- scenes[idx]
      results_var_name <- paste0("Scenario_Medium", scene, "_Fast_Results")
      params_var_name <- paste0("Scenario_Medium", scene, "_Fast_Params")
      
      # Construct the file name dynamically for the results and parameters
      file_name_results <- paste0("~/Transmission_dynamic_model/Output/Scenario_Medium", scene, "_Fast_Results.Rdata")
      file_name_params <- paste0("~/Transmission_dynamic_model/Output/Scenario_Medium", scene, "_Fast_Params.Rdata")
      
      # Save the results and parameters separately
      save(list = "results", file = file_name_results)
      save(list = "params", file = file_name_params)
    }
  }

}

stopCluster(cl) 



