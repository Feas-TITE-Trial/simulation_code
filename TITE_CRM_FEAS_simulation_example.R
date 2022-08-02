# Late onset code for dissertation

# TITE CRM with feasibility

library(dfcrm)
library(poisson)

directory_root <- "~/Desktop/"

source(paste(directory_root,"TITE CRM FEAS CODE/curve_generation_functions.R", sep=""))
source(paste(directory_root,"TITE CRM FEAS CODE/HFD_function.R", sep=""))
source(paste(directory_root,"TITE CRM FEAS CODE/TITE_CRM_FEAS_estimation_function.R", sep=""))
source(paste(directory_root,"TITE CRM FEAS CODE/TITE_CRM_FEAS_scenario_sim.R", sep=""))
source(paste(directory_root,"TITE CRM FEAS CODE/TITE_CRM_FEAS_sim_one_trial.R", sep=""))
source(paste(directory_root,"TITE CRM FEAS CODE/TITE_CRM_stopping_rules.R", sep=""))





# to turn feasibility off set feas_truth to c(1,1,1,1)

scenario_sim1 <- large_simulation_trial_TITE_feas(tox_truth=c(.1,.2,.3,.4),feas_truth=c(.95,.85,.75,.65),skeleton=c(.08,.25,.38,.45),tox_targ=.25,feas_targ=.80,
                                                  pu_tox=.95,pu_feas=.95,max_treat=24,max_extract=30,formal_bayesian=TRUE,model="empiric",intcpt=3,
                                                  prior_feas=c(1,.99,.90,.80,.70),enrollment_dist="poisson",arrival_rate=5,observation_window=10,time_tox_dist="uniform",
                                                  tox_stop=TRUE,partial_obs=FALSE,tox_model_stop="empiric",s2=1.34,init_dose=2,ntrial=2000)
