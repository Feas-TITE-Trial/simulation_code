# Toxicity stopping rules for TITE with feasibility

# partial_obs set to true uses all the as whereas partial obs set to FALSE using only fully observed data


tox_stop_function <- function(skeleton,y_long,dose_long,weight,pu_tox,tox_targ,partial_obs,tox_model_stop,s2){
  
  # skeleton  : crm skeleton
  # y_long    : vector containing toxicitiy indicator for patient i
  # dose_long : vector containing which dose level the ith patient was treated at
  # weight    : vector containing the weights for the i patients so far
  # pu_tox    : upper probability cutoff for toxicity stopping rule
  # tox_targ  : target DLT rate
  # partial obs : logical - TRUE uses partial non DLTs; FALSE uses only fully observed outcomes
  # tox_model_stop: option to change which CRM model the stopping rule is based off of: choices are
  #                 "empiric", "hyperbolic_tangent" , "comp_log_log" , "prod_of_beta" , "probit"
  # s2 : prior variance of parameter
  
  
  
  
  
  weight<- weight[!is.na(weight)]; 
  dose_long<- dose_long[!is.na(dose_long)]; 
  y_long<- y_long[!is.na(y_long)]; 
  
  
  ntreat <- length(dose_long)
  if(partial_obs ==F){
    # y_long <- y_long[weight==1]
    # dose_long <- dose_long[weight==1]
    # weight <- weight[weight==1]
    
    ntreat <- sum(weight==1)
  }
  # print(ntreat)
  # print(weight)
  # print(y_long)
  # print(dose_long)
  if(ntreat>=1){
    
  
  
  likelihood_times_prior <- function(b,skeleton,y_long,dose_long,weight,ntreat,tox_targ,tox_model_stop,s2){
    prior <- exp((-.5*b^2)/s2)
    lik <- 1
    
    
    if(tox_model_stop=="empiric"){
      
    for(i in 1:ntreat){
      x <- skeleton[dose_long[i]] 
      Fx <- x^exp(b)
      lik <- lik*(weight[i]*Fx)^y_long[i] *(1-weight[i]*Fx)^(1-y_long[i])
    }
    
    
    } else if(tox_model_stop=="hyperbolic_tangent"){
      for(i in 1:ntreat){
        x <- skeleton[dose_long[i]] 
        Fx <- ((tanh(x)+1)/2)^exp(b)
        lik <- lik*(weight[i]*Fx)^y_long[i] *(1-weight[i]*Fx)^(1-y_long[i])
      }
    } else if(tox_model_stop=="comp_log_log"){
      for(i in 1:ntreat){
        x <- skeleton[dose_long[i]] 
        Fx <- 1-(1-x)^(exp(b))
        lik <- lik*(weight[i]*Fx)^y_long[i] *(1-weight[i]*Fx)^(1-y_long[i])
      }
    } else if(tox_model_stop=="prod_of_beta"){
      for(i in 1:ntreat){
        x <- skeleton[dose_long[i]] 
        Fx <-1-exp(-b)*(1-x)
        lik <- lik*(weight[i]*Fx)^y_long[i] *(1-weight[i]*Fx)^(1-y_long[i])
      }
    } else if(tox_model_stop=="probit"){
      for(i in 1:ntreat){
        x <- skeleton[dose_long[i]] 
        Fx <-pnorm(b+qnorm(x))
        lik <- lik*(weight[i]*Fx)^y_long[i] *(1-weight[i]*Fx)^(1-y_long[i])
      }
  
    }
    prod <- prior*lik
    return(prod)
  }
  
  if(tox_model_stop == "prod_of_beta"){lower_denom<- log(1-skeleton[1])} else{lower_denom<- -Inf}
    
    
  den <- integrate(likelihood_times_prior,lower=lower_denom,upper=Inf,skeleton,
                   y_long,dose_long,weight,ntreat,tox_targ,tox_model_stop,s2)$value
  
  
  if(tox_model_stop=="empiric"){
  num <- integrate(likelihood_times_prior,lower=-Inf,upper=log(log(tox_targ)/log(skeleton[1])),skeleton,
                   y_long,dose_long,weight,ntreat,tox_targ,tox_model_stop,s2)$value
  } else if (tox_model_stop=="hyperbolic_tangent"){
    num <- integrate(likelihood_times_prior,lower=-Inf,upper=log(log(tox_targ)/((tanh(skeleton[1])+1)/2)),skeleton,
                     y_long,dose_long,weight,ntreat,tox_targ,tox_model_stop,s2)$value
  }  else if (tox_model_stop=="comp_log_log"){
    num <- integrate(likelihood_times_prior,lower=log(log(1-tox_targ)/log(1-skeleton[1])),upper=Inf,skeleton,
                     y_long,dose_long,weight,ntreat,tox_targ,tox_model_stop,s2)$value
  } else if (tox_model_stop=="prod_of_beta"){
    num <- integrate(likelihood_times_prior,lower=log((1-skeleton[1])/(1-tox_targ)),upper=Inf,skeleton,
                     y_long,dose_long,weight,ntreat,tox_targ,tox_model_stop,s2)$value
  } else if (tox_model_stop=="probit"){
    num <- integrate(likelihood_times_prior,lower=qnorm(tox_targ)-qnorm(skeleton[1]),upper=Inf,skeleton,
                     y_long,dose_long,weight,ntreat,tox_targ,tox_model_stop,s2)$value
  }
  
  val <- num/den 
  
  
  if(val>pu_tox){
    dose_1_too_toxic <- TRUE
  } else dose_1_too_toxic <- FALSE
  
  } else if(ntreat<1){val<- 0; dose_1_too_toxic <- FALSE;}
  
  
  
  return(list(dose_1_too_toxic=dose_1_too_toxic,val=val))
  
  
}


