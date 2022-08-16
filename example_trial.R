
# Example Trial Simulation

# Install and load in required packages
library(dfcrm)
library(poisson)
library(nimble)

# Set directory for users computer
directory_root <- "~/Desktop/TITE FEAS/"


# Run required functions to simulate trial
source(paste(directory_root,"hfd_function.R", sep=""))
source(paste(directory_root,"crm_estimation_function.R", sep=""))
source(paste(directory_root,"scenario_sim_function.R", sep=""))
source(paste(directory_root,"trial_sim_function.R", sep=""))
source(paste(directory_root,"stopping_rules.R", sep=""))
source(paste(directory_root,"TITE adaptive weight function.R", sep=""))



# Reproduces the proposed design for scenario 1 in manuscript

results_scenario1_TITE <- large_simulation_trial_TITE_feas(tox_truth=c(.22,.32,.38,.46),feas_truth=c(.89,.70,.37,.25),skeleton=c(.13,.25,.41,.59),tox_targ=.25,feas_targ=.80,
                                                      pu_tox=.90,pu_feas=.90,max_treat=24,max_extract=30,formal_bayesian=TRUE,model="probit_fixed_int",intcpt=3,
                                                      prior_feas=c(1,.80,.60,.40,.20),enrollment_dist="poisson",arrival_rate=5,observation_window=10,time_tox_dist="uniform",
                                                      tox_stop=TRUE,partial_obs=FALSE,tox_model_stop="probit",s2=.74^2,init_dose=1,
                                                      weight_func="linear",allow_dose_skip=F,aplus=1,ntrial=2000)



results_scenario1_TITE$avg_treat
results_scenario1_TITE$fmtd_table_perc






