###############################################################################
# Evaluating the Impact and Cost-Effectiveness of Typhoid Conjugate Vaccine Schedules across Diverse Settings: a Multi-Model Comparison 
# 
# Function to set up model runs based on specified scenarios (incidence, vaccine schedule) under slow-waning assumptions
#
# Created: Catherine Wenger Jan 2025
###############################################################################

run_simulations_parallel <- function(p,scene,arch){
  
  library(RSpectra)
  library(deSolve)
  library(tidyverse)
  library(here)
  
  ## Source relevant scripts -----------------------------------------------------
  source(here("typhoid_ode_slow_waning.R"))                              # The ODE
  
  #########################################################
  # prepare parameters and matrices needed to run the ODE #
  #########################################################
  
  # Demographic variables:
  # B - birth rate I assume
  # mu - death rate
  # N0 - starting population
  # Age groups (al=9)
  #    1 = 0-<9 months
  #    2 = 9-<14 months
  #    3 = 15mo-<2 years
  #    4 = 2-<5 years
  #    5 = 5-<9 years
  #    6 = 10-<15 years
  #    7 = 15-<30 years
  #    8 = 30-<50 years
  #    9 = 50-<70 years
  agepar <- list(B = 0.000438, #per week
                 mu = c(7.027859e-05,  9.418454e-05,  4.345000e-04,  
                        4.117768e-04 ,-2.398326e-04 ,-3.312705e-04, 
                        6.951404e-04  ,3.958299e-04  ,1.329275e-04),
                 N0 = 100000)
  
  
  # Set calibrated parameters based on the incidence setting currently being run
    # b1 - relative risk < 5 yrs old
    # b2 - relative risk for 5-9 yrs old
    # psymp1 - Probability of symptoms (i.e., first infection)
  if(arch == "Medium"){
    archetype = 1                           # set Medium incidence index
    unkp = c(0.2347062,0.56871,0.0283287)   # b1 / b2 / psymp1
  }else if(arch == "High"){
    archetype = 2                           # set High incidence index
    unkp = c(0.2548045,0.9655815,0.1121312) # b1 / b2 / psymp1
  }else if(arch == "VeryHigh"){
    archetype = 3                           # set Very High incidence index
    unkp = c(0.2526783,1.001078,0.6054534)  # b1 / b2 / psymp1
  }  
  
  # Modeling Scenario Set-up ----------------------------------------------------------
  t0=5200             # Length of burn-in period
  tmod = (52*20) + 1  # Model for 20 years post vaccine introduction
  
  tspan=t0+tmod       # Total time span
  al = 9              # The age groups
  
  # Vaccine parameters -----------------------------------------------------------
  ## This specifies the vax coverage for each vaccination scenario
  vaxnon = rbind(matrix(0,tspan, al)) # for no vax, coverage = 0
  ## National campaign vaccination coverage --------------------------------------
  vcn = search_grid[p,'vcn']          # get the stochastic parameter sample for campaign coverage from the search grid for simulation # p
  vcn_r = vcn^(1/4)                   # weekly increase rate (reach maximum coverage by the 4th week after introduction)
  
  ## Routine vaccination coverage ------------------------------------------------
  vcr9 = search_grid[p,'vcr9']    # get the stochastic parameter sample for routine vax @ 9mo from the search grid for simulation # p
  vcr15 = search_grid[p,'vcr15']  # get the stochastic parameter sample for routine vax @ 15mo from the search grid for simulation # p
  vcr2 = search_grid[p,'vcr2']    # get the stochastic parameter sample for routine vax @ 2yr from the search grid for simulation # p
  vcr5 = search_grid[p,'vcr5']    # get the stochastic parameter sample for routine vax @ 5yr from the search grid for simulation # p

  ## Booster vaccination coverage ------------------------------------------------
  vcrb = search_grid[p,'vcrb']    # get the stochastic parameter sample for booster vax @ 5/10yr from the search grid for simulation # p

  ## Vaccine efficacy ------------------------------------------------------------
  VE0 = search_grid[p,'VE0']    # get the stochastic parameter sample for VE @ 0-<2yr from the search grid for simulation # p
  VE1 = search_grid[p,'VE1']    # get the stochastic parameter sample for VE @ 2-<5yr from the search grid for simulation # p
  VE2 = search_grid[p,'VE2']    # get the stochastic parameter sample for VE @ 5+yr from the search grid for simulation # p
  
  ## Rate of waning immunity from vaccination ------------------------------------
  omegav0 = search_grid[p,'omegav0']    # get the stochastic parameter sample for duration of vax protection @ 0-<2yr from the search grid for simulation # p
  omegav1 = search_grid[p,'omegav1']    # get the stochastic parameter sample for duration of vax protection @ 2-<5yr from the search grid for simulation # p
  omegav2 = search_grid[p,'omegav2']    # get the stochastic parameter sample for duration of vax protection @ 5+yr from the search grid for simulation # p
  
  age = c(0.4167,	0.2083,	1.583,	3.5,	7.5,	12.5,	22.5,	40,	60)               # midpoint of each age group
  agelength = c(0.833333333,	0.42,	0.82,	2.93,	5,	5,	15,	20,	20)             # size of age group
  u = c(1/43.3,	1/21.7,	1/42,	1/153,	1/260,	1/260,	1/780,	1/1040,	1/1040)   # aging rate
  ages = c('0-1y','2-5y','5-10y','10-14y','15+y',"All")                         # labels for collapsed age bands (primarily used for fitting purposes)
  
  al = length(age)                                                                  # number of age groups
  agep = c(0.01833,	0.0092,	0.0165,	0.0640,	0.1050,	0.1020,	0.2406,	0.2222,	0.2222) # proportion of population in each age group
  
  N = as.numeric(agepar[3])*agep                                                    # population size in each age group
  B = as.numeric(agepar[1])*cbind(matrix(1,tspan,1) , matrix(0,tspan,al-1))         # Matrix for applying birthrate (B=0 for al>=2)
  N0 = N                                                                            # age-specific starting population size
  mu = as.numeric(agepar[[2]])                                                      # age-specific all-cause mortality rate
  
  omega = 1/104   # rate of waning natural immunity (1/weeks)
  alpha = .001    # proportion of typhoid infections that die
  delta = .25     # rate of recovery from infection (1/weeks)
  
  # Age-group specific: update here for updated age groups
  theta = c(0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.015, 0.066, 0.090)  # proportion of primary typhoid infections who become carriers
  theta2 = rep(0,al)                                                        # proportion of secondary typhoid infections who become carriers

  b1 = unkp[1]              # relative risk of infection for <5 yr olds compared to older ages
  b2 = unkp[2]              # relative risk of infection for 5-9 yr olds compared to older ages
  psymp1 = unkp[3]          # probability of symptoms (first infection)
  epsilon = 0               # rate of returning to fully susceptible state
  r = 1                     # relative infectiousness of non-primary infections
  rC = search_grid[p,'rC']  # get the stochastic parameter sample for relative infectiousness of carriers compared to primary infections from the search grid for simulation # p
  
  # Calculate R0 based on the sampled value of rC and the incidence setting using fit equations
  if (archetype == 1){
    R0 = 2.986 + 18.299*rC + -10.847*rC^2 # Medium
  }else if(archetype == 2){
    R0 = 2.341 + 34.351*rC + -18.846*rC^2 # High
  }else if(archetype == 3){
    R0 = 2.484 + 45.756*rC + -20.057*rC^2 # Very High
  }
  
  ###################################################
  # Create transmission matrix using fit parameters #
  ###################################################
  
  # b1 applies to age groups 1 & 2 & 3, so they are collapsed into 1 age group 
  # in the following 4 lines in preparation of calculating the transmission matrix
  alfit = 7                                                           # number of age groups for the fit parameters
  agepfit = c(0.044, 0.0640,	0.1050,	0.1020,	0.2406,	0.2222,	0.2222) # proportion of population in each age group (Equal for each group)
  mufit = c(((mu[1]*10/24)+(mu[2]*5/24)+(mu[3]*9/24)),mu[4:9])        # age-specific mortality rate
  thetafit = c(.003, .003, .003, .003, .015, .066, .090)              # proportion of primary typhoid infections who become carriers
  
  # person-to-person transmission rate calculation 
  betapfit = R0*(delta+mean(mufit))/
    as.numeric(1 + rC*delta*(agepfit%*%(thetafit/mean(mufit))))*
    rbind(b1*rep(1,alfit),b2*rep(1,alfit), matrix(1,alfit-2,alfit))*alfit/
    max(eigen(t( rbind(b1*rep(1,alfit),b2*rep(1,alfit), matrix(1,alfit-2,alfit))))$values)
  
  # Age-group specific: update here to expand the matrix to match the 9 age groups
  betap = rbind(betapfit[1,],betapfit[1,], betapfit)
  betap = cbind(betap[,1],betap[,1], betap)
  
  # starting population matrix to provide for the ODE
  y0 = c(round(.89*N0)-10,    # [1] S1        - susceptible
         10*rep(1,al),        # [2] I1        - infected
         round(0*N0),         # [3] R         - recovered
         round(.1*N0)-10,     # [4] S2        - susceptible (prior infection)
         10*rep(1,al),        # [5] I2        - infected (prior infection)
         round(.01*N0),       # [6] C         - carrier
         round(0*N0),         # [7] SV1       - susceptible, vaccination
         round(0*N0),         # [8] IV1       - infected, vaccination
         round(0*N0),         # [9] RV        - recovered, vaccination
         round(0*N0),         # [10] SV2      - susceptible, vaccination (prior infection)
         round(0*N0),         # [11] IV2      - infected, vaccination (prior infection)
         round(0*N0),         # [12] C2       - carrier, vaccinated (prior infection)
         round(0*N0),         # [13] V1 <2    - vaccinated at age 0-<2y, protected
         round(0*N0),         # [14] V2 <2    - vaccinated at age 0-<2y, protected (prior infection)
         round(0*N0),         # [15] V1 2-5   - vaccinated at age 2-<5y, protected
         round(0*N0),         # [16] V2 2-5   - vaccinated at age 2-<5y, protected (prior infection)
         round(0*N0),         # [17] V1 >= 5  - vaccinated at age >5y, protected
         round(0*N0),         # [18] V2 >= 5  - vaccinated at age >5y, protected (prior infection)
         round(0*N0),         # [19] Cumulative: Campaign doses
         round(0*N0),         # [20] Cumulative: Boosters doses
         round(0*N0),         # [21] Cumulative: Routine doses
         10*rep(1,al),        # [22] Cumulative: I1 primary infections
         round(0*N0))         # [23] Cumulative: IV1 primary infections (vaccinated)
  
  # combine set parameters to be given to the ODE solver
  pars <- list(al=al,
               epsilon=epsilon,
               alpha= alpha,
               delta = delta,
               omega=omega,
               gamma=gamma,
               r=r,
               rC=rC,
               betap=betap,
               B=B,
               u=u,
               mu=mu,
               theta=theta,
               theta2=  theta2,
               omegav0 = omegav0,
               omegav1 = omegav1,
               omegav2 = omegav2,
               VE0 = VE0,
               VE1=VE1,
               VE2=VE2,
               END_OF_BURN_IN=t0)
  
  #####################################################################################
  # prepare vaccination matrix based on strategy being implemented in this simulation #
  #####################################################################################
  
  if(scene == "A"){
    #   A - No vaccination
    pars <- c(pars,
              list(vaxr = vaxnon,   # routine coverage is zero
                   vaxb = vaxnon,   # booster coverage is zero
                   vaxnc = vaxnon)) # campaign coverage is zero
  }else if(scene == "B"){
    #   B - Routine @ 9mo without catch-up
    vaxr9 = rbind(matrix(0, t0, al),                                       # Pre-vax era
                  matrix(c(vcr9,0,0,0,0,0,0,0,0), tmod, al, byrow = TRUE)) # vaccination as individuals age into the second age group
    pars <- c(pars,
              list(vaxr = vaxr9,    # routine coverage is set for 9 months
                   vaxb = vaxnon,   # booster coverage is zero
                   vaxnc = vaxnon)) # campaign coverage is zero
  }else if(scene == "C"){
    #   C - Routine @ 9mo with catch-up
    vaxr9 = rbind(matrix(0, t0, al),                                                                # Pre-vax era
                  matrix(c(vcr9,0,0,0,0,0,0,0,0), tmod, al, byrow = TRUE))                          # vaccination as individuals age into the second age group
    vaxnc9 = rbind(matrix(0, t0, al),                                                               # Pre-vax era
                   matrix(c(0,vcn_r^4,vcn_r^4,vcn_r^4,vcn_r^4,vcn_r^4,0,0,0), 1, al, byrow = TRUE), # week 1 campaign coverage
                   matrix(c(0,vcn_r^3,vcn_r^3,vcn_r^3,vcn_r^3,vcn_r^3,0,0,0), 1, al, byrow = TRUE), # week 2 campaign coverage
                   matrix(c(0,vcn_r^2,vcn_r^2,vcn_r^2,vcn_r^2,vcn_r^2,0,0,0), 1, al, byrow = TRUE), # week 3 campaign coverage
                   matrix(c(0,vcn_r,vcn_r,vcn_r,vcn_r,vcn_r,0,0,0), 1, al, byrow = TRUE),           # week 4 campaign coverage (full coverage)
                   matrix(c(0,0,0,0,0,0,0,0,0), tmod-4, al, byrow = TRUE))                          # campaign coverage ended
    pars <- c(pars,
              list(vaxr = vaxr9,    # routine coverage is set for 9 months
                   vaxb = vaxnon,   # booster coverage is zero
                   vaxnc = vaxnc9)) # campaign coverage is set for 9 months - 15 years
  } else if(scene == "D"){
    #   D - Routine @ 15mo without catch-up
    vaxr15 = rbind(matrix(0, t0, al),                                        # Pre-vax era
                   matrix(c(0,vcr15,0,0,0,0,0,0,0), tmod, al, byrow = TRUE)) # vaccination as individuals age into the third age group
    pars <- c(pars,
              list(vaxr = vaxr15,   # routine coverage is set for 15 months
                   vaxb = vaxnon,   # booster coverage is zero
                   vaxnc = vaxnon)) # campaign coverage is zero
  } else if(scene == "E"){
    #   E - Routine @ 15mo with catch-up
    vaxr15 = rbind(matrix(0, t0, al),                                                          # Pre-vax era
                   matrix(c(0,vcr15,0,0,0,0,0,0,0), tmod, al, byrow = TRUE))                   # vaccination as individuals age into the third age group 
    vaxnc15 = rbind(matrix(0, t0, al),                                                         # Pre-vax era
                    matrix(c(0,0,vcn_r^4,vcn_r^4,vcn_r^4,vcn_r^4,0,0,0), 1, al, byrow = TRUE), # week 1 campaign coverage
                    matrix(c(0,0,vcn_r^3,vcn_r^3,vcn_r^3,vcn_r^3,0,0,0), 1, al, byrow = TRUE), # week 2 campaign coverage
                    matrix(c(0,0,vcn_r^2,vcn_r^2,vcn_r^2,vcn_r^2,0,0,0), 1, al, byrow = TRUE), # week 3 campaign coverage
                    matrix(c(0,0,vcn_r,vcn_r,vcn_r,vcn_r,0,0,0), 1, al, byrow = TRUE),         # week 4 campaign coverage (full coverage)
                    matrix(c(0,0,0,0,0,0,0,0,0), tmod-4, al, byrow = TRUE))                    # campaign coverage ended
    pars <- c(pars,
              list(vaxr = vaxr15,    # routine coverage is set for 15 months
                   vaxb = vaxnon,    # booster coverage is zero
                   vaxnc = vaxnc15)) # campaign coverage is set for 15 months - 15 years
  } else if(scene == "F"){
    #   F - Routine @ 2yrs without catch-up
    vaxr2 = rbind(matrix(0, t0, al),                                       # Pre-vax era
                  matrix(c(0,0,vcr2,0,0,0,0,0,0), tmod, al, byrow = TRUE)) # vaccination as individuals age into the fourth age group
    pars <- c(pars,
              list(vaxr = vaxr2,    # routine coverage is set for 2 years
                   vaxb = vaxnon,   # booster coverage is zero
                   vaxnc = vaxnon)) # campaign coverage is zero
  } else if(scene == "G"){
    #   G - Routine @ 2yrs with catch-up
    vaxr2 = rbind(matrix(0, t0, al),                                                    # Pre-vax era
                  matrix(c(0,0,vcr2,0,0,0,0,0,0), tmod, al, byrow = TRUE))              # vaccination as individuals age into the fourth age group
    vaxnc2 = rbind(matrix(0, t0, al),                                                   # Pre-vax era
                   matrix(c(0,0,0,vcn_r^4,vcn_r^4,vcn_r^4,0,0,0), 1, al, byrow = TRUE), # week 1 campaign coverage
                   matrix(c(0,0,0,vcn_r^3,vcn_r^3,vcn_r^3,0,0,0), 1, al, byrow = TRUE), # week 2 campaign coverage
                   matrix(c(0,0,0,vcn_r^2,vcn_r^2,vcn_r^2,0,0,0), 1, al, byrow = TRUE), # week 3 campaign coverage
                   matrix(c(0,0,0,vcn_r,vcn_r,vcn_r,0,0,0), 1, al, byrow = TRUE),       # week 4 campaign coverage (full coverage)
                   matrix(c(0,0,0,0,0,0,0,0,0), tmod-4, al, byrow = TRUE))              # campaign coverage ended
    pars <- c(pars,
              list(vaxr = vaxr2,    # routine coverage is set for 2 years
                   vaxb = vaxnon,   # booster coverage is zero
                   vaxnc = vaxnc2)) # campaign coverage is set for 2 - 15 years
  } else if(scene == "H"){
    #   H - Routine @ 5yrs without catch-up
    vaxr5 = rbind(matrix(0, t0, al),                                        # Pre-vax era
                  matrix(c(0,0,0,vcr5,0,0,0,0,0), tmod, al, byrow = TRUE))  # vaccination as individuals age into the fifth age group
    pars <- c(pars,
              list(vaxr = vaxr5,    # routine coverage is set for 5 years
                   vaxb = vaxnon,   # booster coverage is zero
                   vaxnc = vaxnon)) # campaign coverage is zero
  } else if(scene == "I"){
    #   I - Routine @ 5yrs with catch-up
    vaxr5 = rbind(matrix(0, t0, al),                                              # Pre-vax era
                  matrix(c(0,0,0,vcr5,0,0,0,0,0), tmod, al, byrow = TRUE))        # vaccination as individuals age into the fifth age group
    vaxnc5 = rbind(matrix(0, t0, al),                                             # Pre-vax era
                   matrix(c(0,0,0,0,vcn_r^4,vcn_r^4,0,0,0), 1, al, byrow = TRUE), # week 1 campaign coverage
                   matrix(c(0,0,0,0,vcn_r^3,vcn_r^3,0,0,0), 1, al, byrow = TRUE), # week 2 campaign coverage
                   matrix(c(0,0,0,0,vcn_r^2,vcn_r^2,0,0,0), 1, al, byrow = TRUE), # week 3 campaign coverage
                   matrix(c(0,0,0,0,vcn_r,vcn_r,0,0,0), 1, al, byrow = TRUE),     # week 4 campaign coverage (full coverage)
                   matrix(c(0,0,0,0,0,0,0,0,0), tmod-4, al, byrow = TRUE))        # campaign coverage ended
    pars <- c(pars,
              list(vaxr = vaxr5,    # routine coverage is set for 5 years
                   vaxb = vaxnon,   # booster coverage is zero
                   vaxnc = vaxnc5)) # campaign coverage is set for 5 - 15 years
  } else if(scene == "J"){
    #   J - Routine @ 9mo without catch-up + 5yrs Booster
    vaxr9 = rbind(matrix(0, t0, al),                                                # Pre-vax era
                  matrix(c(vcr9,0,0,0,0,0,0,0,0), tmod, al, byrow = TRUE))          # vaccination as individuals age into the second age group
    vaxb5 = rbind(matrix(0, t0 + (52*4), al),                                       # Pre-vax era
                  matrix(c(0,0,0,vcrb,0,0,0,0,0), tmod - (52*4), al, byrow = TRUE)) # vaccination as individuals age into the fifth age group
    pars <- c(pars,
              list(vaxr = vaxr9,    # routine coverage is set for 9 months
                   vaxb = vaxb5,    # booster coverage is set for 5 years
                   vaxnc = vaxnon)) # campaign coverage is zero
  } else if(scene == "K"){
    #   K - Routine @ 9mo with catch-up + 5yrs Booster
    vaxr9 = rbind(matrix(0, t0, al),                                                                # Pre-vax era
                  matrix(c(vcr9,0,0,0,0,0,0,0,0), tmod, al, byrow = TRUE))                          # vaccination as individuals age into the second age group
    vaxb5 = rbind(matrix(0, t0 + (52*4), al),                                                       # Pre-vax era
                  matrix(c(0,0,0,vcrb,0,0,0,0,0), tmod - (52*4), al, byrow = TRUE))                 # vaccination as individuals age into the fifth age group
    vaxnc9 = rbind(matrix(0, t0, al),                                                               # Pre-vax era
                   matrix(c(0,vcn_r^4,vcn_r^4,vcn_r^4,vcn_r^4,vcn_r^4,0,0,0), 1, al, byrow = TRUE), # week 1 campaign coverage
                   matrix(c(0,vcn_r^3,vcn_r^3,vcn_r^3,vcn_r^3,vcn_r^3,0,0,0), 1, al, byrow = TRUE), # week 2 campaign coverage
                   matrix(c(0,vcn_r^2,vcn_r^2,vcn_r^2,vcn_r^2,vcn_r^2,0,0,0), 1, al, byrow = TRUE), # week 3 campaign coverage
                   matrix(c(0,vcn_r,vcn_r,vcn_r,vcn_r,vcn_r,0,0,0), 1, al, byrow = TRUE),           # week 4 campaign coverage (full coverage)
                   matrix(c(0,0,0,0,0,0,0,0,0), tmod-4, al, byrow = TRUE))                          # campaign coverage ended
    pars <- c(pars,
              list(vaxr = vaxr9,    # routine coverage is set for 9 months
                   vaxb = vaxb5,    # booster coverage is set for 5 years
                   vaxnc = vaxnc9)) # campaign coverage is set for 9 months - 15 years
  } else if(scene == "L"){
    #   L - Routine @ 15mo without catch-up + 5yrs Booster
    vaxr15 = rbind(matrix(0, t0, al),                                                # Pre-vax era
                   matrix(c(0,vcr15,0,0,0,0,0,0,0), tmod, al, byrow = TRUE))         # vaccination as individuals age into the third age group 
    vaxb5 = rbind(matrix(0, t0 + (52*4), al),                                        # Pre-vax era
                  matrix(c(0,0,0,vcrb,0,0,0,0,0), tmod - (52*4), al, byrow = TRUE))  # vaccination as individuals age into the fifth age group
    pars <- c(pars,
              list(vaxr = vaxr15,   # routine coverage is set for 15 months
                   vaxb = vaxb5,    # booster coverage is set for 5 years
                   vaxnc = vaxnon)) # campaign coverage is zero
  } else if(scene == "M"){
    #   M - Routine @ 15mo with catch-up + 5yrs Booster
    vaxr15 = rbind(matrix(0, t0, al),                                                          # Pre-vax era
                   matrix(c(0,vcr15,0,0,0,0,0,0,0), tmod, al, byrow = TRUE))                   # vaccination as individuals age into the third age group 
    vaxb5 = rbind(matrix(0, t0 + (52*4), al),                                                  # Pre-vax era
                  matrix(c(0,0,0,vcrb,0,0,0,0,0), tmod - (52*4), al, byrow = TRUE))            # vaccination as individuals age into the fifth age group
    vaxnc15 = rbind(matrix(0, t0, al),                                                         # Pre-vax era
                    matrix(c(0,0,vcn_r^4,vcn_r^4,vcn_r^4,vcn_r^4,0,0,0), 1, al, byrow = TRUE), # week 1 campaign coverage
                    matrix(c(0,0,vcn_r^3,vcn_r^3,vcn_r^3,vcn_r^3,0,0,0), 1, al, byrow = TRUE), # week 2 campaign coverage
                    matrix(c(0,0,vcn_r^2,vcn_r^2,vcn_r^2,vcn_r^2,0,0,0), 1, al, byrow = TRUE), # week 3 campaign coverage
                    matrix(c(0,0,vcn_r,vcn_r,vcn_r,vcn_r,0,0,0), 1, al, byrow = TRUE),         # week 4 campaign coverage (full coverage)
                    matrix(c(0,0,0,0,0,0,0,0,0), tmod-4, al, byrow = TRUE))                    # campaign coverage ended
    pars <- c(pars,
              list(vaxr = vaxr15,    # routine coverage is set for 15 months
                   vaxb = vaxb5,     # booster coverage is set for 5 years
                   vaxnc = vaxnc15)) # campaign coverage is set for 15 months - 15 years
  } else if(scene == "N"){
    #   N - Routine @ 9mo without catch-up + 5 & 10yrs Booster
    vaxr9 = rbind(matrix(0, t0, al),                                                      # Pre-vax era
                  matrix(c(vcr9,0,0,0,0,0,0,0,0), tmod, al, byrow = TRUE))                # vaccination as individuals age into the second age group 
    vaxb5_10 = rbind(matrix(0, t0 + (52*4), al),                                          # Pre-vax era
                     matrix(c(0,0,0,vcrb,0,0,0,0,0), (52*5), al, byrow = TRUE),           # vaccination as individuals age into the fifth age group
                     matrix(c(0,0,0,vcrb,vcrb,0,0,0,0), tmod - (52*9), al, byrow = TRUE)) # vaccination as individuals age into the fifth age group
    pars <- c(pars,
              list(vaxr = vaxr9,    # routine coverage is set for 9 months
                   vaxb = vaxb5_10, # booster coverage is set for 5&10 years
                   vaxnc = vaxnon)) # campaign coverage is zero
  } else if(scene == "O"){
    #   O - Routine @ 9mo with catch-up + 5 & 10yrs Booster
    vaxr9 = rbind(matrix(0, t0, al),                                                                # Pre-vax era
                  matrix(c(vcr9,0,0,0,0,0,0,0,0), tmod, al, byrow = TRUE))                          # vaccination as individuals age into the second age group 
    vaxb5_10 = rbind(matrix(0, t0 + (52*4), al),                                                    # Pre-vax era
                     matrix(c(0,0,0,vcrb,0,0,0,0,0), (52*5), al, byrow = TRUE),                     # vaccination as individuals age into the fifth age group
                     matrix(c(0,0,0,vcrb,vcrb,0,0,0,0), tmod - (52*9), al, byrow = TRUE))           # vaccination as individuals age into the fifth age group
    vaxnc9 = rbind(matrix(0, t0, al),                                                               # Pre-vax era
                   matrix(c(0,vcn_r^4,vcn_r^4,vcn_r^4,vcn_r^4,vcn_r^4,0,0,0), 1, al, byrow = TRUE), # week 1 campaign coverage
                   matrix(c(0,vcn_r^3,vcn_r^3,vcn_r^3,vcn_r^3,vcn_r^3,0,0,0), 1, al, byrow = TRUE), # week 2 campaign coverage
                   matrix(c(0,vcn_r^2,vcn_r^2,vcn_r^2,vcn_r^2,vcn_r^2,0,0,0), 1, al, byrow = TRUE), # week 3 campaign coverage
                   matrix(c(0,vcn_r,vcn_r,vcn_r,vcn_r,vcn_r,0,0,0), 1, al, byrow = TRUE),           # week 4 campaign coverage (full coverage)
                   matrix(c(0,0,0,0,0,0,0,0,0), tmod-4, al, byrow = TRUE))                          # campaign coverage ended
    pars <- c(pars,
              list(vaxr = vaxr9,    # routine coverage is set for 9 months
                   vaxb = vaxb5_10, # booster coverage is set for 5&10 years
                   vaxnc = vaxnc9)) # campaign coverage is set for 9 months - 15 years
  } else if(scene == "P"){
    #   P - Routine @ 15mo without catch-up + 5 & 10yrs Booster
    vaxr15 = rbind(matrix(0, t0, al),                                                     # Pre-vax era
                   matrix(c(0,vcr15,0,0,0,0,0,0,0), tmod, al, byrow = TRUE))              # vaccination as individuals age into the third age group  
    vaxb5_10 = rbind(matrix(0, t0 + (52*4), al),                                          # Pre-vax era
                     matrix(c(0,0,0,vcrb,0,0,0,0,0), (52*5), al, byrow = TRUE),           # vaccination as individuals age into the fifth age group
                     matrix(c(0,0,0,vcrb,vcrb,0,0,0,0), tmod - (52*9), al, byrow = TRUE)) # vaccination as individuals age into the fifth age group
    pars <- c(pars,
              list(vaxr = vaxr15,   # routine coverage is set for 15 months
                   vaxb = vaxb5_10, # booster coverage is set for 5&10 years
                   vaxnc = vaxnon)) # campaign coverage is zero
  } else if(scene == "Q"){
    #   Q - Routine @ 15mo with catch-up + 5 & 10yrs Booster
    vaxr15 = rbind(matrix(0, t0, al),                                                          # Pre-vax era
                   matrix(c(0,vcr15,0,0,0,0,0,0,0), tmod, al, byrow = TRUE))                   # vaccination as individuals age into the third age group  
    vaxb5_10 = rbind(matrix(0, t0 + (52*4), al),                                               # Pre-vax era
                     matrix(c(0,0,0,vcrb,0,0,0,0,0), (52*5), al, byrow = TRUE),                # vaccination as individuals age into the fifth age group
                     matrix(c(0,0,0,vcrb,vcrb,0,0,0,0), tmod - (52*9), al, byrow = TRUE))      # vaccination as individuals age into the fifth age group
    vaxnc15 = rbind(matrix(0, t0, al),                                                         # Pre-vax era
                    matrix(c(0,0,vcn_r^4,vcn_r^4,vcn_r^4,vcn_r^4,0,0,0), 1, al, byrow = TRUE), # week 1 campaign coverage
                    matrix(c(0,0,vcn_r^3,vcn_r^3,vcn_r^3,vcn_r^3,0,0,0), 1, al, byrow = TRUE), # week 2 campaign coverage
                    matrix(c(0,0,vcn_r^2,vcn_r^2,vcn_r^2,vcn_r^2,0,0,0), 1, al, byrow = TRUE), # week 3 campaign coverage
                    matrix(c(0,0,vcn_r,vcn_r,vcn_r,vcn_r,0,0,0), 1, al, byrow = TRUE),         # week 4 campaign coverage (full coverage)
                    matrix(c(0,0,0,0,0,0,0,0,0), tmod-4, al, byrow = TRUE))                    # campaign coverage ended
    pars <- c(pars,
              list(vaxr = vaxr15,    # routine coverage is set for 15 months
                   vaxb = vaxb5_10,  # booster coverage is set for 5&10 years
                   vaxnc = vaxnc15)) # campaign coverage is set for 15 months - 15 years
  }

  # run the ODE using the previously set parameters across the defined time span
  results <- as.data.frame(ode(y = y0,                 # starting population for each age group/compartment
                               t=1:tspan,              # defined time span
                               parms = pars,           # all set parameters
                               func = typhoid_ode_slow, # ODE function to be called by the solver
                               method = "lsoda",       # solving method
                               rtol = 1e-6))           # tolerance
  
  # Return a list with the simulation results AND the stochastic (non-fixed) parameter values
  return(list(
    #results = results,
    results = results[(99*52):(tspan-1),],  # the main simulation results
    non_fixed_parms = data.frame(VE0 = VE0,
                                 VE1 = VE1,
                                 VE2 = VE2,
                                 omegav0 = omegav0,
                                 omegav1 = omegav1,
                                 omegav2 = omegav2,
                                 vcn=vcn,
                                 vcr9=vcr9,
                                 vcr15=vcr15,
                                 vcr2=vcr2,
                                 vcr5=vcr5,
                                 vcrb=vcrb,
                                 rC = rC)
  ))
  
  gc()
}
