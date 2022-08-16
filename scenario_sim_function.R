# TITE CRM FEAS simulation for an entire scenario

large_simulation_trial_TITE_feas <- function(tox_truth,feas_truth,skeleton,tox_targ,feas_targ,pu_tox,pu_feas,max_treat,max_extract,formal_bayesian,model,intcpt,
                                             prior_feas,enrollment_dist,arrival_rate,observation_window,time_tox_dist,tox_stop,partial_obs,
                                             tox_model_stop,s2,init_dose,weight_func,allow_dose_skip,wshape,wscale,aplus,restrict,ntrial){
  ndose <- length(tox_truth)
  nmat <- ymat <- matrix(NA,nrow=ntrial,ncol=ndose)
  mtd_vec <- fmtd_vec <- stop_feas_vec <- stop_tox_vec <- duration_vec <- HFD_final_vec <-  rep(NA,ntrial)
  
  
  all_res <- NULL
  
  for(i in 1:ntrial){
    
    res <- simulate_trial_TITE_feas(tox_truth,feas_truth,skeleton,tox_targ,feas_targ,pu_tox,pu_feas,max_treat,max_extract,formal_bayesian,model,intcpt,
                                    prior_feas,enrollment_dist,arrival_rate,observation_window,time_tox_dist,tox_stop,partial_obs,tox_model_stop,s2,init_dose,weight_func,
                                    allow_dose_skip,wshape,wscale,aplus,restrict)
    
    
    
    
    all_res[[i]] <- res
    nmat[i,] <- res$n
    ymat[i,] <- res$y
    mtd_vec[i] <- res$crm_final_est
    fmtd_vec[i] <- res$FMTD_est
    stop_feas_vec[i] <- res$stop_feas
    stop_tox_vec[i] <- res$stop_tox
    duration_vec[i] <- res$duration
    HFD_final_vec[i] <- res$HFD_final
    #print(i)
    
  }
  avg_treat <- round(apply(nmat,MARGIN=2,FUN=mean),3)
  avg_tox <- round(apply(ymat,MARGIN=2,FUN=mean),3)
  
  mean_duration <- mean(duration_vec)
  
  mtd_table_perc <- as.numeric(table(factor(mtd_vec,levels = 0:ndose)))/ntrial*100; names(mtd_table_perc) <- as.character(0:ndose);
  fmtd_table_perc <- as.numeric(table(factor(fmtd_vec,levels = 0:ndose)))/ntrial*100; names(fmtd_table_perc) <- as.character(0:ndose);
  hfd_table_perc <- as.numeric(table(factor(HFD_final_vec,levels = 0:ndose)))/ntrial*100; names(hfd_table_perc) <- as.character(0:ndose);
  perc_stop_feas <- sum(stop_feas_vec)/ntrial*100
  perc_stop_tox <- sum(stop_tox_vec)/ntrial*100
  
  return(list(avg_treat=avg_treat,avg_tox=avg_tox,
              mtd_vec=mtd_vec,fmtd_vec=fmtd_vec,stop_feas_vec=stop_feas_vec,stop_tox_vec=stop_tox_vec,
              duration_vec=duration_vec,nmat=nmat,ymat=ymat,HFD_final_vec=HFD_final_vec,mean_duration=mean_duration,
              mtd_table_perc=mtd_table_perc,fmtd_table_perc=fmtd_table_perc,hfd_table_perc=hfd_table_perc,
              perc_stop_feas=perc_stop_feas,perc_stop_tox=perc_stop_tox))
  
}

