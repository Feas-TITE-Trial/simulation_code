# this function simulates a single trial



simulate_trial_TITE_feas <- function(tox_truth,feas_truth,skeleton,tox_targ,feas_targ,pu_tox,pu_feas,max_treat,max_extract,formal_bayesian,model,intcpt,
                                     prior_feas,enrollment_dist,arrival_rate,observation_window,time_tox_dist,tox_stop,partial_obs,tox_model_stop,s2,init_dose){
  
  ndose <- length(tox_truth)  # Number of dose levels
  dose_levels <- 1:ndose      # contains the possible doses
  n <- y <- rep(0,ndose)      # Count and outcomes at each dose level
  tox_risk <- runif(n=max_treat) # Patient toxicity risks for patients who are treatable at some dose level
  feas_prob <- c(1-feas_truth[1], abs(diff(feas_truth)), feas_truth[ndose]) # probability vector for highest feasible doses
  patient_feas <- sample(x=c(0:ndose),size = max_extract,prob = feas_prob,replace = T) # patients highest feasible dose levels
  
  ### for testing tite crm ####
  #patient_feas <- rep(ndose,max_extract)
  
  patient_ID <- 1:max_treat # Identifies patients who are put on trial
  X <- rep(0,ndose+1)
  crm_estimate <- rep(NA,max_treat)
  
  weight <- dose_long <- y_long <- feas_ind <- tox_ind <- event_time <- treat_time <- rep(NA,max_treat);
  
  ### Enrollment times and time to toxicity###
  # Default options
  if(missing(enrollment_dist)) enrollment_dist <- "poisson"
  if(missing(arrival_rate)) arrival_rate <- 5
  if(missing(observation_window)) observation_window <- 10
  if(missing(time_tox_dist)) time_tox_dist <- "uniform"
  if(missing(tox_stop)) tox_stop <- F
  if(missing(partial_obs)) partial_obs <- T
  if(missing(tox_model_stop)) tox_model_stop <- "empiric"
  
  if(enrollment_dist == "poisson"){
    enrollment_time <- hpp.event.times(rate=arrival_rate/observation_window,num.events = max_extract)
  }
  if(time_tox_dist == "uniform") {tox_time <- enrollment_time + runif(n=max_extract,min=0,max=observation_window)
  } else if(time_tox_dist=="weibull"){
    w_times <- rweibull(n=max_extract,shape=5,scale=5)
    w_times[w_times>observation_window] <- observation_window-.01
    tox_time <- enrollment_time + w_times
  }
  
  stop_tox <- stop_feas <- F # indicators for whether the trial should be stopped
  
  init_dose <- init_dose  # initial dose specified in function argument
  
  curr_target_dose <- init_dose # this is the current targeted dose
  
  i <- 0
  patient_ID <- 0
  
  while(i <max_extract){
    
    # Enroll Patient and determine their HFD
    i <- i + 1
    hfd_curr <- patient_feas[i]
    X[hfd_curr+1] <- X[hfd_curr+1] + 1
    if(hfd_curr == 0){ # check that dose level 1 is still feasible
      global_HFD <- HFD_check(X,feas_targ,pu_feas,prior_feas)
      if(global_HFD==0){
        stop_feas <- T        # Feasibility stopping rule
        break
      } else{ # move on to next patient if current patient not feasible but globally dose 1 still feasible
        next
      }
    } else{ # current patient is feasible at some dose level
      
      ###### update current patient ouctomes ######
      
      
      curr_time <- enrollment_time[i]
      weight <- curr_time - treat_time
      y_long[curr_time>event_time & tox_ind == 1] <- 1
      
      w <- weight/observation_window
      w[w>1] <- 1
      
      #######
      
      #### check stopping rule for toxicity
      if(patient_ID > 1){
        if(tox_stop_function(skeleton=skeleton,y_long=y_long,dose_long=dose_long,
                             weight=w,pu_tox=pu_tox,tox_targ=tox_targ,partial_obs,tox_model_stop=tox_model_stop,s2=s2)$dose_1_too_toxic){
          stop_tox <- T
          break
        }
      }
      
      ####
      
      # if we don't stop the trial due to excessive toxicity then we enroll a patient
      patient_ID <- patient_ID+1
      
      
      #### update mtd estimate according to the TITE CRM with current data
      if(patient_ID>1){
        est_curr <- get_tox_estimate(skeleton = skeleton,ylong=y_long,doselong = dose_long,
                                     weight =w,ntreat=sum(n),target = tox_targ,formal_bayesian=formal_bayesian,model=model,s2=s2)
        crm_estimate[i] <- est_curr$model_mtd
        curr_target_dose <- crm_estimate[i] 
      }
      
      
      ### treat current patient with updated target dose
      
      dose_long[patient_ID] <- min(curr_target_dose,hfd_curr)
      n[dose_long[patient_ID]] <- n[dose_long[patient_ID]] + 1
      y_long[patient_ID] <- 0
      #weight[patient_ID] <- 0
      treat_time[patient_ID] <- curr_time
      tox_ind[patient_ID] <- as.numeric(tox_risk[patient_ID] < tox_truth[dose_long[patient_ID]]) # this determines if a patient will have a DLT
      y[dose_long[patient_ID]] <- y[dose_long[patient_ID]] + tox_ind[patient_ID]
      if(tox_ind[patient_ID]==1){event_time[patient_ID] <- tox_time[i]} else event_time[patient_ID] <- treat_time[patient_ID] + observation_window
      if(sum(n)>=max_treat){break} #stop trial when it is fully enrolled
    }
    
    
    
    
    
  }
  # compute final mtd at the conclusion of the trial based upon all the data 
  y_long <- tox_ind
  if(stop_feas==F  & stop_tox == F){
    est_final <- get_tox_estimate(skeleton = skeleton,ylong=tox_ind,doselong = dose_long,
                                  weight =rep(1,max_treat),ntreat=sum(n),target = tox_targ,
                                  formal_bayesian=formal_bayesian,model = model,s2=s2)
    crm_final_est <- est_final$model_mtd
    HFD_final <- HFD_check(X,feas_targ,pu_feas,prior_feas)
    FMTD_est <- min(crm_final_est,HFD_final)
  } else{
    crm_final_est <- FMTD_est <- HFD_final <- 0
  }
  
  if(sum(n)>0) {duration <- max(event_time,na.rm=T)} else duration<- 0
  
  
  
  
  return(list(n=n,y=y,dose_long=dose_long,y_long=y_long,treat_time=treat_time,event_time=event_time,crm_final_est=crm_final_est,FMTD_est=FMTD_est,
              stop_tox=stop_tox,stop_feas=stop_feas,duration=duration,HFD_final=HFD_final))
  
}