# TITE adaptive weight function from Cheung and Chappell

weight_calc_function <- function(obs_window,DLT_status,time_on_trial){
  
  time_on_trial <-time_on_trial[!is.na(time_on_trial)]
  DLT_status <-DLT_status[!is.na(DLT_status)]
  time_on_trial[time_on_trial>obs_window] <- obs_window
  
  npatients <- length(time_on_trial)
  
  ordered_times <- sort(c(0,time_on_trial[as.logical(DLT_status)],obs_window))

  weight_vec <- rep(NA, npatients)
  
  z <- sum(DLT_status)
  
  

  
  
  
  for(i in 1:npatients){
    
    if(DLT_status[i]==1 | time_on_trial[i]>=obs_window){ 
      weight_vec[i] <- 1
    } else{
      u <- time_on_trial[i]
      if(z==0){
         
        k <- 0
        
      } else{
        ugreater <- u >= ordered_times
        k <- sum(ugreater)-1
      }
      weight_vec[i] <- (k/(z+1)) + (1/(z+1))*((u-ordered_times[k+1])/(ordered_times[k+2]-ordered_times[k+1]))
      
    }
      
   
    
    
    
  }
  
  return(list(weight_vec=weight_vec))
}





