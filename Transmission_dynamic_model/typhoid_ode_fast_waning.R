###############################################################################
# Evaluating the Impact and Cost-Effectiveness of Typhoid Conjugate Vaccine Schedules across Diverse Settings: a Multi-Model Comparison 
# 
# Function containing the Ordinary Differential Equations used in the ODE solver in the function "parallelization_fast"
#     There are two compartments for each vaccination class the fast-waning scenario follows a gamma distribution (shape = 2)
#     First compartment represents the early phase of waning, the second represents a later phase of waning
#
# Created: Catherine Wenger Jan 2025
###############################################################################

typhoid_ode_fast <- function(t,y,parameters){
  
  al = parameters$al                          # number of age groups
  
  # Set up states
  dydt = vector(len = 29*al)

  # Parameters (alphabetical)
  alpha = parameters$alpha                    # proportion of typhoid infections that die
  B = parameters$B                            # Birth rate
  betap = parameters$betap                    # Transmission rate matrix
  delta = parameters$delta                    # Rate of recovery from typhoid infection
  epsilon = parameters$epsilon                # Rate of returning to a fully susceptible state (0)
  END_OF_BURN_IN = parameters$END_OF_BURN_IN  # Timepoint at the end of the burn-in period
  mu = parameters$mu                          # age-specific all-cause death rate
  omega = parameters$omega                    # duration of narutal immunity
  omegav0_1 = parameters$omegav0_1            # Duration of vaccine protection for individuals vaccinated at 0-<2 years (early phase of waning)
  omegav0_2 = parameters$omegav0_2            # Duration of vaccine protection for individuals vaccinated at 0-<2 years (later phase of waning)
  omegav1_1 = parameters$omegav1_1            # Duration of vaccine protection for individuals vaccinated at 2-<5 years (early phase of waning)
  omegav1_2 = parameters$omegav1_2            # Duration of vaccine protection for individuals vaccinated at 2-<5 years (later phase of waning)
  omegav2_1 = parameters$omegav2_1            # Duration of vaccine protection for individuals vaccinated at 5+ years (early phase of waning)
  omegav2_2 = parameters$omegav2_2            # Duration of vaccine protection for individuals vaccinated at 5+ years (later phase of waning)
  r = parameters$r                            # relative infectiousness of non-primary infections (1, they are equal)
  rC = parameters$rC                          # relatie infectiousness of chronic carriers compared to primary infections
  theta = parameters$theta                    # proportion of primary typhoid infections who become carriers
  theta2 = parameters$theta2                  # proportion of secondary typhoid infections who become carriers (0)
  u = parameters$u                            # Aging rate (dependent on length of time in each age band)
  vaxr = parameters$vaxr                      # Coverage for routine vaccination
  vaxnc = parameters$vaxnc                    # Coverage for vaccine campaigns
  vaxb = parameters$vaxb                      # Coverage for booster vaccination
  VE0 = parameters$VE0                        # VE for individuals vaccinated at 0-<2 years
  VE1 = parameters$VE1                        # VE for individuals vaccinated at 2-<5 years
  VE2 = parameters$VE2                        # VE for individuals vaccinated at 5+ years
  t <- min(max(1, round(t)), tspan)           # Update the time step of the simulation from the ODE solver

  # Get the sum of entire population for time step t
  N = sum(y[1:(length(y)-(5*al))])
  
  # Get indicies for compartments contributing to the force of infection calculation
  I1_indices = (al + 1):(2*al)  
  I2_indices = (4*al + 1):(5*al)
  C_indices = (5*al + 1):(6*al)
  IV1_indices = (7*al + 1):(8*al)
  IV2_indices = (10*al + 1):(11*al)
  CV_indices = (11*al + 1):(12*al)
  # Force of infection for short-term cycle
  #I1 + r*I2 + rC*C + IV1 + r*IV2 + rC*CV / sum of population
  lambdap = betap %*% ((y[I1_indices] + r*y[I2_indices] + rC*y[C_indices] + y[IV1_indices] + r*y[IV2_indices] + rC*y[CV_indices])/N)
  
  lambdaw = 0 # no water transmission
  
  # >>> GENERAL STRUCTURE OF THE ODES <<<
  ## [1] Disease-related dynamics
  ## [2] Vaccination through campaigns
  ## [3] Aging processess + routine and booster vaccinations (which are part of aging)

  # Loop through each age groups (al)
  for(i in 1:al){

        # Initialize each state holder, so it is easier to keep track in the ODE equations
        S1 = y[i]                 # susceptible
        I1 = y[i+al]              # infected
        R = y[i+2*al]             # recovered
        S2 = y[i+3*al]            # susceptible (prior infection)
        I2 = y[i+4*al]            # infected (prior infection)
        C = y[i+5*al]             # carrier
        SV1 = y[i+6*al]           # vaccinated, susceptible
        IV1 = y[i+7*al]           # vaccinated, infected
        RV = y[i+8*al]            # vaccinated, recovered
        SV2 = y[i+9*al]           # vaccinated, susceptible (prior infection)
        IV2 = y[i+10*al]          # vaccinated, infected (prior infection)
        CV = y[i+11*al]           # vaccinated, carrier
        V_UNDER_2_V1_1 = y[i+12*al] # includes routine vaccinees and campaign vaccinees under 2 at time of campaign (early phase of waning)
        V_UNDER_2_V1_2 = y[i+13*al] # includes routine vaccinees and campaign vaccinees under 2 at time of campaign (later phase of waning)
        V_UNDER_2_V2_1 = y[i+14*al] # includes routine vaccinees and campaign vaccinees under 2 at time of campaign (prior infection, early phase of waning)
        V_UNDER_2_V2_2 = y[i+15*al] # includes routine vaccinees and campaign vaccinees under 2 at time of campaign (prior infection, later phase of waning)
        V_2TO5_V1_1 = y[i+16*al]    # includes routine vaccinees and campaign vaccinees 2-5 at time of campaign (early phase of waning)
        V_2TO5_V1_2 = y[i+17*al]    # includes routine vaccinees and campaign vaccinees 2-5 at time of campaign (later phase of waning)
        V_2TO5_V2_1 = y[i+18*al]    # includes routine vaccinees and campaign vaccinees 2-5 at time of campaign (prior infection, early phase of waning)
        V_2TO5_V2_2 = y[i+19*al]    # includes routine vaccinees and campaign vaccinees 2-5 at time of campaign (prior infection, later phase of waning)
        V_OVER_5_V1_1 = y[i+20*al]  # includes routine/booster recipients and campaign recipients over 5 at time of campaign (early phase of waning)
        V_OVER_5_V1_2 = y[i+21*al]  # includes routine/booster recipients and campaign recipients over 5 at time of campaign (later phase of waning)
        V_OVER_5_V2_1 = y[i+22*al]  # includes routine/booster recipients and campaign recipients over 5 at time of campaign (prior infection, early phase of waning)
        V_OVER_5_V2_2 = y[i+23*al]  # includes routine/booster recipients and campaign recipients over 5 at time of campaign (prior infection, later phase of waning)

        # >>> BLOCK [1] DISEASE-TRANSMISSION
        # >>> UNVACCINATED FOLKS <<<
        
        dS1 = B[round(t),i]*N + epsilon*S2 - (lambdap[i]+lambdaw)*S1 - (u[i]+mu[i])*S1
        dI1 = (lambdap[i]+lambdaw)*S1 - delta*(1-alpha-theta[i])*I1 - delta*theta[i]*I1 - delta*alpha*I1 - (u[i]+mu[i])*I1
        dR = delta*(1-alpha-theta[i])*I1 + delta*(1-theta2[i])*I2 - omega*R - (u[i]+mu[i])*R
        dS2 = omega*R - epsilon*S2 - (lambdap[i]+lambdaw)*S2 - (u[i]+mu[i])*S2
        dI2 = (lambdap[i]+lambdaw)*S2 - delta*I2 - (u[i]+mu[i])*I2
        dC =  delta*(theta[i]*I1 + theta2[i]*I2) - (u[i]+mu[i])*C

        # >>> VACCINATED FOLKS <<<
        # Disease classes among previously vaccinated

        dSV1 = omegav0_1*V_UNDER_2_V1_1 + omegav0_2*V_UNDER_2_V1_2 + 
                        omegav1_1*V_2TO5_V1_1 + omegav1_2*V_2TO5_V1_2 + 
                        omegav2_1*V_OVER_5_V1_1 + omegav2_2*V_OVER_5_V1_2 + 
                        epsilon*SV2 - (lambdap[i]+lambdaw)*SV1 - (u[i]+mu[i])*SV1
        dIV1 = (lambdap[i]+lambdaw)*SV1 + 
                        (lambdap[i]+lambdaw)*(1-VE0)*V_UNDER_2_V1_1 + (lambdap[i]+lambdaw)*(1-VE0)*V_UNDER_2_V1_2 + 
                        (lambdap[i]+lambdaw)*(1-VE1)*V_2TO5_V1_1 + (lambdap[i]+lambdaw)*(1-VE1)*V_2TO5_V1_2 +
                        (lambdap[i]+lambdaw)*(1-VE2)*V_OVER_5_V1_1 + (lambdap[i]+lambdaw)*(1-VE2)*V_OVER_5_V1_2 + 
                        - delta*(1-alpha-theta[i])*IV1 - delta*theta[i]*IV1 - delta*alpha*IV1 - (u[i]+mu[i])*IV1
        dRV = delta*(1-alpha-theta[i])*IV1 + delta*(1-theta2[i])*IV2 - omega*RV - (u[i]+mu[i])*RV
        dSV2 = omegav0_1*V_UNDER_2_V2_1 + omegav0_2*V_UNDER_2_V2_2 + 
                        omegav1_1*V_2TO5_V2_1 + omegav1_2*V_2TO5_V2_2 + 
                        omegav2_1*V_OVER_5_V2_1 + omegav2_2*V_OVER_5_V2_2 + 
                        omega*RV - epsilon*SV2 - (lambdap[i]+lambdaw)*SV2 - (u[i]+mu[i])*SV2
        dIV2 = (lambdap[i]+lambdaw)*SV2 + 
                        (lambdap[i]+lambdaw)*(1-VE0)*V_UNDER_2_V2_1 + (lambdap[i]+lambdaw)*(1-VE0)*V_UNDER_2_V2_2 + 
                        (lambdap[i]+lambdaw)*(1-VE1)*V_2TO5_V2_1 + (lambdap[i]+lambdaw)*(1-VE1)*V_2TO5_V2_2 + 
                        (lambdap[i]+lambdaw)*(1-VE2)*V_OVER_5_V2_1 + (lambdap[i]+lambdaw)*(1-VE2)*V_OVER_5_V2_2 + 
                        - delta*IV2 - (u[i]+mu[i])*IV2
        dCV = delta*(theta[i]*IV1 + theta2[i]*IV2) - (u[i]+mu[i])*CV

        # VACCINATION CLASSESS
        # there are two compartments for each vaccination class the fast-waning scenario follows a gamma distribution (shape = 2)
        # "_1" designation indicates early phase of waning, "_2" designation indicates later phase of waning
        # If vaccinated at a time when they were <2 years of age
        dV_UNDER_2_V1_1 = -omegav0_1*V_UNDER_2_V1_1 - (u[i]+mu[i])*V_UNDER_2_V1_1 - (lambdap[i]+lambdaw)*(1-VE0)*V_UNDER_2_V1_1 
        dV_UNDER_2_V1_2 = -omegav0_2*V_UNDER_2_V1_2 - (u[i]+mu[i])*V_UNDER_2_V1_2 - (lambdap[i]+lambdaw)*(1-VE0)*V_UNDER_2_V1_2 
        
        dV_UNDER_2_V2_1 = -omegav0_1*V_UNDER_2_V2_1 - (u[i]+mu[i])*V_UNDER_2_V2_1 - (lambdap[i]+lambdaw)*(1-VE0)*V_UNDER_2_V2_1
        dV_UNDER_2_V2_2 = -omegav0_2*V_UNDER_2_V2_2 - (u[i]+mu[i])*V_UNDER_2_V2_2 - (lambdap[i]+lambdaw)*(1-VE0)*V_UNDER_2_V2_2
        
        # If vaccinated at a time when they were 2-<5 years of age
        dV_2TO5_V1_1 = -omegav1_1*V_2TO5_V1_1 - (u[i]+mu[i])*V_2TO5_V1_1 - (lambdap[i]+lambdaw)*(1-VE1)*V_2TO5_V1_1 
        dV_2TO5_V1_2 = -omegav1_2*V_2TO5_V1_2 - (u[i]+mu[i])*V_2TO5_V1_2 - (lambdap[i]+lambdaw)*(1-VE1)*V_2TO5_V1_2 
        
        dV_2TO5_V2_1 = -omegav1_1*V_2TO5_V2_1 - (u[i]+mu[i])*V_2TO5_V2_1 - (lambdap[i]+lambdaw)*(1-VE1)*V_2TO5_V2_1
        dV_2TO5_V2_2 = -omegav1_2*V_2TO5_V2_2 - (u[i]+mu[i])*V_2TO5_V2_2 - (lambdap[i]+lambdaw)*(1-VE1)*V_2TO5_V2_2

        # If vaccinated at a time when they were 5 or older
        dV_OVER_5_V1_1 = -omegav2_1*V_OVER_5_V1_1 - (u[i]+mu[i])*V_OVER_5_V1_1 - (lambdap[i]+lambdaw)*(1-VE2)*V_OVER_5_V1_1
        dV_OVER_5_V1_2 = -omegav2_2*V_OVER_5_V1_2 - (u[i]+mu[i])*V_OVER_5_V1_2 - (lambdap[i]+lambdaw)*(1-VE2)*V_OVER_5_V1_2
        
        dV_OVER_5_V2_1 = -omegav2_1*V_OVER_5_V2_1 - (u[i]+mu[i])*V_OVER_5_V2_1 - (lambdap[i]+lambdaw)*(1-VE2)*V_OVER_5_V2_1
        dV_OVER_5_V2_2 = -omegav2_2*V_OVER_5_V2_2 - (u[i]+mu[i])*V_OVER_5_V2_2 - (lambdap[i]+lambdaw)*(1-VE2)*V_OVER_5_V2_2

        # >>> COUNTERS for doses and symptomatic infections <<<
        dCuml_Vax_Camp_Doses = 0
        dCuml_Vax_Booster_Doses = 0
        dCuml_Vax_Routine_Doses = 0
        dCuml_I1 = (lambdap[i] + lambdaw) * S1
        dCuml_IV1 = (lambdap[i] + lambdaw) * SV1 + 
                        (lambdap[i]+lambdaw)*(1-VE0)*V_UNDER_2_V1_1 + (lambdap[i]+lambdaw)*(1-VE0)*V_UNDER_2_V1_2 + 
                        (lambdap[i]+lambdaw)*(1-VE1)*V_2TO5_V1_1 + (lambdap[i]+lambdaw)*(1-VE1)*V_2TO5_V1_2 + 
                        (lambdap[i]+lambdaw)*(1-VE2)*V_OVER_5_V1_1 + (lambdap[i]+lambdaw)*(1-VE2)*V_OVER_5_V1_2 


        # >>> BLOCK [2] VACCINATION: CAMPAIGNS
        # If out of the burn-in period and into the vaccination time
        # This when we begin vaccination
        # leaky vaccine protection: all vaccinated individuals move to vaccine classes, but have a chance of protection in those classes related to VE
        
        if(round(t) > END_OF_BURN_IN){
            # individuals removed from their unvaccinated disease compartments based on vax coverage
            dS1 = dS1 - S1*(vaxnc[round(t),i])
            dI1 = dI1 - I1*(vaxnc[round(t),i])
            dR = dR - R*(vaxnc[round(t),i])
            dS2 = dS2 - S2*(vaxnc[round(t),i])
            dI2 = dI2 - I2*(vaxnc[round(t),i])
            dC = dC - C*(vaxnc[round(t),i])
            # infected individuals and carriers move to vaccinated, unprotected compartments 
            dIV1 = dIV1 + I1*(vaxnc[round(t),i])
            dIV2= dIV2 + I2*(vaxnc[round(t),i])
            dCV = dCV + C*(vaxnc[round(t),i])

            ## ! PEOPLE MOVING INTO VACCINATION WILL DEPEND ON THE AGE AT WHICH THEY ARE VACCCINATED
            ## If i (age group counter) is 1, 2, or 3 (0-2 yrs) then they move into <2 yr group
            #  If i is 3 (i.e.2-<4y) they move into the 2-5 yr group
            #  Else if i is 4-al (i.e.y) they move into the >5 yr group
            #  vaccinated individuals only enter "_1" designated vaccine compartments
            if(i %in% c(1,2,3)){
              dSV1 = dSV1
              dRV = dRV
              dSV2 = dSV2

              # Changes in V1 <2 yr group
              # people who were vaccinated in the catch-up campaign moving into vaccinated class based on prior infection
              dV_UNDER_2_V1_1 = dV_UNDER_2_V1_1 + S1*(vaxnc[round(t),i])      # no prior infection 
              dV_UNDER_2_V2_1 = dV_UNDER_2_V2_1 + (S2+R)*(vaxnc[round(t),i])  # prior infection 
                  

            } else if(i %in% c(4)){
              dSV1 = dSV1
              dRV = dRV
              dSV2 = dSV2
              
              # Changes in V1 2-5 yr group
              # people who were vaccinated in the catch-up campaign moving into vaccinated class based on prior infection
              dV_2TO5_V1_1 = dV_2TO5_V1_1 + S1*(vaxnc[round(t),i])      # no prior infection 
              dV_2TO5_V2_1 = dV_2TO5_V2_1 + (S2+R)*(vaxnc[round(t),i])  # prior infection
              
            } else if(i %in% c(5:al)){
              dSV1 = dSV1
              dRV = dRV
              dSV2 = dSV2

              # Changes in V1 >5 yr group
              # people who were vaccinated in the catch-up campaign moving into vaccinated class based on prior infection
              dV_OVER_5_V1_1 = dV_OVER_5_V1_1 + S1*(vaxnc[round(t),i])      # no prior infection
              dV_OVER_5_V2_1 = dV_OVER_5_V2_1 + (S2+R)*(vaxnc[round(t),i])  # prior infection
            }

            # Cumulative vaccinations through campaign doses
            dCuml_Vax_Camp_Doses = dCuml_Vax_Camp_Doses + ((S1 + I1 + R + S2 + I2 + C) * vaxnc[round(t),i])

          }

        # >>> BLOCK [3] AGING + BOOSTER OR ROUTINE VACCINES
        # NOTE: Since we dealt with AGING OUT in the first block, this only handles where people go (i.e., AGING IN)
        # THEREFORE, we only need to start beyond the first age group
            if (i>1){
              # size of each compartment at the previous time step
              S1_prior = y[i-1]
              I1_prior = y[i+al-1]
              R_prior = y[i+2*al-1]
              S2_prior = y[i+3*al-1]
              I2_prior = y[i+4*al-1]
              C_prior = y[i+5*al-1]
              SV1_prior = y[i+6*al-1]
              IV1_prior = y[i+7*al-1]
              RV_prior = y[i+8*al-1]
              SV2_prior = y[i+9*al-1]
              IV2_prior = y[i+10*al-1]
              CV_prior = y[i+11*al-1]
              V_UNDER_2_V1_prior_1 = y[i+12*al-1]
              V_UNDER_2_V1_prior_2 = y[i+13*al-1]
              V_UNDER_2_V2_prior_1 = y[i+14*al-1]
              V_UNDER_2_V2_prior_2 = y[i+15*al-1]
              V_2TO5_V1_prior_1 = y[i+16*al-1]
              V_2TO5_V1_prior_2 = y[i+17*al-1]
              V_2TO5_V2_prior_1 = y[i+18*al-1]
              V_2TO5_V2_prior_2 = y[i+19*al-1]
              V_OVER_5_V1_prior_1 = y[i+20*al-1]
              V_OVER_5_V1_prior_2 = y[i+21*al-1]
              V_OVER_5_V2_prior_1 = y[i+22*al-1]
              V_OVER_5_V2_prior_2 = y[i+23*al-1]

              # >>> | ROUTINE DOSES | <<<
              # >>> UNVACCINATED FOLKS <<<
              ## This is the folks who don't get vaccinated at all =  1 - Vaccination Coverage (routine)
              dS1 = dS1 + u[i-1]*S1_prior*(1-(vaxr[round(t),(i-1)]))
              dI1 = dI1 + u[i-1]*I1_prior*(1-(vaxr[round(t),(i-1)]))
              dR = dR + u[i-1]*R_prior*(1-(vaxr[round(t),(i-1)]))
              dS2 = dS2 + u[i-1]*S2_prior*(1-(vaxr[round(t),(i-1)]))
              dI2 = dI2 + u[i-1]*I2_prior*(1-(vaxr[round(t),(i-1)]))
              dC = dC + u[i-1]*C_prior*(1-(vaxr[round(t),(i-1)]))

              # >>> VACCINATED FOLKS <<<
              ## This is compromised of folks who get vaccinated, but infecteds/carriers do not gain immunity from the vaccine
              dSV1 = dSV1
              dIV1 = dIV1 + u[i-1]*I1_prior*vaxr[round(t),(i-1)] 
              dRV = dRV
              dSV2 = dSV2
              dIV2 = dIV2 + u[i-1]*I2_prior*vaxr[round(t),(i-1)] 
              dCV = dCV + u[i-1]*C_prior*vaxr[round(t),(i-1)] 

              # Aging of vaccinated but not protected groups
              dIV1 = dIV1 + u[i-1]*IV1_prior
              dIV2 = dIV2 + u[i-1]*IV2_prior
              dCV = dCV + u[i-1]*CV_prior

              # Aging of vaccine-protected groups (out of younger vaccine-protected compartments, into older vaccine-protected compartments)
              dV_UNDER_2_V1_1 = dV_UNDER_2_V1_1 + u[i-1]*V_UNDER_2_V1_prior_1
              dV_UNDER_2_V1_2 = dV_UNDER_2_V1_2 + u[i-1]*V_UNDER_2_V1_prior_2
              
              dV_UNDER_2_V2_1 = dV_UNDER_2_V2_1 + u[i-1]*V_UNDER_2_V2_prior_1
              dV_UNDER_2_V2_2 = dV_UNDER_2_V2_2 + u[i-1]*V_UNDER_2_V2_prior_2
              
              dV_2TO5_V1_1 = dV_2TO5_V1_1 + u[i-1]*V_2TO5_V1_prior_1
              dV_2TO5_V1_2 = dV_2TO5_V1_2 + u[i-1]*V_2TO5_V1_prior_2
              
              dV_2TO5_V2_1 = dV_2TO5_V2_1 + u[i-1]*V_2TO5_V2_prior_1
              dV_2TO5_V2_2 = dV_2TO5_V2_2 + u[i-1]*V_2TO5_V2_prior_2
              
              dV_OVER_5_V1_1 = dV_OVER_5_V1_1 + u[i-1]*V_OVER_5_V1_prior_1
              dV_OVER_5_V1_2 = dV_OVER_5_V1_2 + u[i-1]*V_OVER_5_V1_prior_2
              
              dV_OVER_5_V2_1 = dV_OVER_5_V2_1 + u[i-1]*V_OVER_5_V2_prior_1
              dV_OVER_5_V2_2 = dV_OVER_5_V2_2 + u[i-1]*V_OVER_5_V2_prior_2

              # vaccinating individuals as they age in to the age group for routine vaccination
              #  vaccinated individuals only enter "_1" designated vaccine compartments
              if(i %in% c(2,3)){
                # Individuals from S/R compartments aging into vaccine eligible age group and getting vaccinated (<2 yr group)
                dV_UNDER_2_V1_1 = dV_UNDER_2_V1_1 + u[i-1]*S1_prior*(vaxr[round(t),(i-1)])
                dV_UNDER_2_V2_1 = dV_UNDER_2_V2_1 + u[i-1]*(S2_prior+R_prior)*(vaxr[round(t),(i-1)])
              } else if(i %in% c(4)){
                # Individuals from S/R compartments aging into vaccine eligible age group and getting vaccinated (2-5 yr group)
                dV_2TO5_V1_1 = dV_2TO5_V1_1 + u[i-1]*S1_prior*(vaxr[round(t),(i-1)])
                dV_2TO5_V2_1 = dV_2TO5_V2_1 + u[i-1]*(S2_prior+R_prior)*(vaxr[round(t),(i-1)])
              } else if(i %in% c(5:al)){
                # Individuals from S/R compartments aging into vaccine eligible age group and getting vaccinated (>5 yr group)
                dV_OVER_5_V1_1 = dV_OVER_5_V1_1 + u[i-1]*S1_prior*(vaxr[round(t),(i-1)])
                dV_OVER_5_V2_1 = dV_OVER_5_V2_1 + u[i-1]*(S2_prior+R_prior)*(vaxr[round(t),(i-1)])
              }

              # Cumulative vaccinations - routine
              dCuml_Vax_Routine_Doses = dCuml_Vax_Routine_Doses + u[i-1]*(S1_prior + I1_prior + R_prior + S2_prior + I2_prior + C_prior)*(vaxr[round(t),(i-1)])

              # >>> | BOOSTER DOSES | <<<
              ## This is for the folks who don't get vaccinated at all =  1 - Vaccination Coverage (booster) 
              dSV1 = dSV1 +  u[i-1]*SV1_prior*(1-(vaxb[round(t),(i-1)]))  # Those who don't receive a booster (based on coverage)
              dRV = dRV +  u[i-1]*RV_prior*(1-(vaxb[round(t),(i-1)]))     # Those who don't receive a booster (based on coverage)
              dSV2 = dSV2 + u[i-1]*SV2_prior*(1-(vaxb[round(t),(i-1)]))   # Those who don't receive a booster (based on coverage)
              
              # vaccinating individuals as they age in to the age group for booster vaccination
              # only individuals who have been previously vaccinated can get a booster dose
              #  vaccinated individuals only enter "_1" designated vaccine compartments
              if(i %in% c(2,3)){
                # Individuals from S/R compartments aging into vaccine eligible age group and getting boosted (<2 yr group)
                dV_UNDER_2_V1_1 = dV_UNDER_2_V1_1 + u[i-1]*SV1_prior*(vaxb[round(t),(i-1)])
                dV_UNDER_2_V2_1 = dV_UNDER_2_V2_1 + u[i-1]*(SV2_prior + RV_prior)*(vaxb[round(t),(i-1)])
                } else if(i %in% c(4)){
                  # Individuals from S/R compartments aging into vaccine eligible age group and getting boosted (2-5 yr group)
                dV_2TO5_V1_1 = dV_2TO5_V1_1 + u[i-1]*SV1_prior*(vaxb[round(t),(i-1)])
                dV_2TO5_V2_1 = dV_2TO5_V2_1 + u[i-1]*(SV2_prior + RV_prior)*(vaxb[round(t),(i-1)])
                } else if(i %in% c(5:al)){
                  # Individuals from S/R compartments aging into vaccine eligible age group and getting boosted (>5 yr group)
                dV_OVER_5_V1_1 = dV_OVER_5_V1_1 + u[i-1]*SV1_prior*(vaxb[round(t),(i-1)])
                dV_OVER_5_V2_1 = dV_OVER_5_V2_1 + u[i-1]*(SV2_prior + RV_prior)*(vaxb[round(t),(i-1)])
                }
              # Cumulative vaccinations booster
                dCuml_Vax_Booster_Doses = dCuml_Vax_Booster_Doses + u[i-1]*(SV1_prior + IV1_prior + RV_prior + SV2_prior + IV2_prior + CV_prior + 
                                            V_UNDER_2_V1_prior_1 + V_UNDER_2_V2_prior_1 + V_2TO5_V1_prior_1 + V_2TO5_V2_prior_1 + V_OVER_5_V1_prior_1 + V_OVER_5_V2_prior_1)*(vaxb[round(t),(i-1)])

            }
        # adding the compartment sizes for the current time step to the overall simulation dataframe
        dydt[i] = dS1
        dydt[i+al] = dI1
        dydt[i+2*al] = dR
        dydt[i+3*al] = dS2
        dydt[i+4*al] = dI2
        dydt[i+5*al] = dC
        dydt[i+6*al] = dSV1
        dydt[i+7*al] = dIV1
        dydt[i+8*al] = dRV
        dydt[i+9*al] = dSV2
        dydt[i+10*al] = dIV2
        dydt[i+11*al] = dCV
        dydt[i+12*al] = dV_UNDER_2_V1_1
        dydt[i+13*al] = dV_UNDER_2_V1_2
        dydt[i+14*al] = dV_UNDER_2_V2_1
        dydt[i+15*al] = dV_UNDER_2_V2_2
        dydt[i+16*al] = dV_2TO5_V1_1
        dydt[i+17*al] = dV_2TO5_V1_2
        dydt[i+18*al] = dV_2TO5_V2_1
        dydt[i+19*al] = dV_2TO5_V2_2
        dydt[i+20*al] = dV_OVER_5_V1_1
        dydt[i+21*al] = dV_OVER_5_V1_2
        dydt[i+22*al] = dV_OVER_5_V2_1
        dydt[i+23*al] = dV_OVER_5_V2_2
        dydt[i+24*al] = dCuml_Vax_Camp_Doses
        dydt[i+25*al] = dCuml_Vax_Booster_Doses
        dydt[i+26*al] = dCuml_Vax_Routine_Doses
        dydt[i+27*al] = dCuml_I1
        dydt[i+28*al] = dCuml_IV1

  }
  # return the compartment data to the solver
  list(dydt)
}
