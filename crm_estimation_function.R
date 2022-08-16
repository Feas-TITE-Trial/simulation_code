# time to event trial - toxicity evaluation 


get_tox_estimate <- function(skeleton,ylong,doselong,weight,ntreat,target,formal_bayesian,model,intcpt,s2){
  
  # Input descriptions
  # skeleton: skeleton for the bayesian crm
  # ylong   : vector containing if the ith patient has experienced a DLT
  # doselong: vector containing which dose level the ith patient was treated at
  # weight  : vector containing the weights for the i patients so far
  # ntreat  : scalar containing number of patients treated so far
  # target  : target toxicity probability
  # formal_bayesian : logical option : TRUE estimates the mean of the model; False uses the plug in estimate (see page 20 of CRM text)
  # model   : option to specify form of the CRM. Choices  are: 
  #           "empiric", "logistic_fixed_int", "logistic_fixed_slope", "hyperbolic_tangent",
  #           "comp_log_log", "prod_of_beta", "probit_fixed_int"," probit_fixed_slope"
  # intcpt  : intercept of the model - defaults to 3 if nothing is entered
  # s2      : prior variance of the parameter

  
  if(missing(model)){model<- "empiric"}
  if(missing(intcpt)){intcpt<- 3}
  #print(ylong); print(doselong); print(weight); print(ntreat);
  
  
  if(model=="empiric"){
    likelihood_times_prior <- function(b,skeleton,ylong,doselong,weight,ntreat,target,s2){
      prior <- exp((-.5*b^2)/s2)
      lik <- 1
      
      for(i in 1:ntreat){
        x <- skeleton[doselong[i]]
        Fx <- x^exp(b)
        lik <- lik*(weight[i]*Fx)^ylong[i] *(1-weight[i]*Fx)^(1-ylong[i])

      }
      prod <- prior*lik
      return(prod)
    }
    
    likelihood_times_prior_mean <- function(b,skeleton,ylong,doselong,weight,ntreat,target,s2,j){
      Fp <- skeleton[j]^exp(b)
      return(Fp*likelihood_times_prior(b,skeleton,ylong,doselong,weight,ntreat,target,s2))
  
      
    }
    

    
    
  } else if(model=="logistic_fixed_int"){
    a0 <- intcpt
    likelihood_times_prior <- function(b,skeleton,ylong,doselong,weight,ntreat,target,s2){
      prior <- exp((-.5*b^2)/s2)
      lik <- 1

      for(i in 1:ntreat){
        x <- skeleton[doselong[i]] 
        logitk <- log(x/(1-x))
        Fx <- (exp(a0+exp(b)*(logitk-a0)))/(1+(exp(a0+exp(b)*(logitk-a0))))
        lik <- lik*(weight[i]*Fx)^ylong[i] *(1-weight[i]*Fx)^(1-ylong[i])
      }
      prod <- prior*lik
      return(prod)
    }
    
    likelihood_times_prior_mean <- function(b,skeleton,ylong,doselong,weight,ntreat,target,s2,j){
      a0 <- intcpt
      logitp <- log(skeleton[j]/(1-skeleton[j]))
      Fp <- (exp(a0+exp(b)*(logitp-a0)))/(1+(exp(a0+exp(b)*(logitp-a0))))
      

      return(Fp*likelihood_times_prior(b,skeleton,ylong,doselong,weight,ntreat,target,s2))
      
    }
    
    } else if(model=="logistic_fixed_slope"){ # this onw isnt working
      likelihood_times_prior <- function(b,skeleton,ylong,doselong,weight,ntreat,target,s2){
        prior <- exp((-.5*b^2)/s2)
        lik <- 1
        
        for(i in 1:ntreat){
          x <- skeleton[doselong[i]] 
          Fx <- (exp(b)*x)/(1-x+exp(b)*x)

          lik <- lik*(weight[i]*Fx)^ylong[i] *(1-weight[i]*Fx)^(1-ylong[i])
        }
        prod <- prior*lik
        return(prod)
      }
      
      likelihood_times_prior_mean <- function(b,skeleton,ylong,doselong,weight,ntreat,target,s2,j){

        Fp  <- (exp(b)*skeleton[j])/(1-skeleton[j]+exp(b)*skeleton[j])
        
        return(Fp*likelihood_times_prior(b,skeleton,ylong,doselong,weight,ntreat,target,s2))
        
      }
      
      
      
  
    

    
  } else if(model=="hyperbolic_tangent"){
    likelihood_times_prior <- function(b,skeleton,ylong,doselong,weight,ntreat,target,s2){
      prior <- exp((-.5*b^2)/s2)
      lik <- 1
      
      for(i in 1:ntreat){
        x <- skeleton[doselong[i]] 
        Fx <- ((tanh(x)+1)/2)^exp(b)
        lik <- lik*(weight[i]*Fx)^ylong[i] *(1-weight[i]*Fx)^(1-ylong[i])
      }
      prod <- prior*lik
      return(prod)
    }
    
    likelihood_times_prior_mean <- function(b,skeleton,ylong,doselong,weight,ntreat,target,s2,j){
      Fp <- ((tanh(skeleton[j])+1)/2)^exp(b)
      
     
      return(Fp*likelihood_times_prior(b,skeleton,ylong,doselong,weight,ntreat,target,s2))
      
    }
    

    
  } else if(model=="comp_log_log"){
    likelihood_times_prior <- function(b,skeleton,ylong,doselong,weight,ntreat,target,s2){
      prior <- exp((-.5*b^2)/s2)
      lik <- 1
      
      for(i in 1:ntreat){
        x <- skeleton[doselong[i]] 
        Fx <- 1-(1-x)^(exp(b))
        lik <- lik*(weight[i]*Fx)^ylong[i] *(1-weight[i]*Fx)^(1-ylong[i])
      }
      prod <- prior*lik
      return(prod)
    }
    
    likelihood_times_prior_mean <- function(b,skeleton,ylong,doselong,weight,ntreat,target,s2,j){
      Fp <- 1-(1-skeleton[j])^(exp(b))
      return(Fp*likelihood_times_prior(b,skeleton,ylong,doselong,weight,ntreat,target,s2))
      
    }
    

    
  } else if(model=="prod_of_beta"){
    likelihood_times_prior <- function(b,skeleton,ylong,doselong,weight,ntreat,target,s2){
      prior <- exp((-.5*b^2)/s2)
      lik <- 1
      
      for(i in 1:ntreat){
        x <- skeleton[doselong[i]] 
        Fx <- 1-(exp(-1*b)*(1-x))
        lik <- lik*(weight[i]*Fx)^ylong[i] *(1-weight[i]*Fx)^(1-ylong[i])
      }
      prod <- prior*lik
      return(prod)
    }
    
    likelihood_times_prior_mean <- function(b,skeleton,ylong,doselong,weight,ntreat,target,s2,j){
      Fp <- 1-exp(-1*b)*(1-skeleton[j])
      return(Fp*likelihood_times_prior(b,skeleton,ylong,doselong,weight,ntreat,target,s2))
      
    }
    

  }  else if(model=="probit_fixed_int"){
    a0 <- intcpt
    likelihood_times_prior <- function(b,skeleton,ylong,doselong,weight,ntreat,target,s2){
      prior <- exp((-.5*b^2)/s2)
      lik <- 1
      
      for(i in 1:ntreat){
        x <- skeleton[doselong[i]] 
        Fx <- pnorm(a0+exp(b)*(qnorm(x)-a0))
        lik <- lik*(weight[i]*Fx)^ylong[i] *(1-weight[i]*Fx)^(1-ylong[i])
      }
      prod <- prior*lik
      return(prod)
    }
    likelihood_times_prior_mean <- function(b,skeleton,ylong,doselong,weight,ntreat,target,s2,j){
      Fp <- pnorm(a0+exp(b)*(qnorm(skeleton[j])-a0))
      return(Fp*likelihood_times_prior(b,skeleton,ylong,doselong,weight,ntreat,target,s2))
      
    }
      
    } else if(model=="probit_fixed_slope"){
   
      likelihood_times_prior <- function(b,skeleton,ylong,doselong,weight,ntreat,target,s2){
        prior <- exp((-.5*b^2)/s2)
        lik <- 1
        
        for(i in 1:ntreat){
          x <- skeleton[doselong[i]] 
          Fx <- pnorm(b+qnorm(x))
          
          lik <- lik*(weight[i]*Fx)^ylong[i] *(1-weight[i]*Fx)^(1-ylong[i])
        }
        prod <- prior*lik
        return(prod)
      }
    
    likelihood_times_prior_mean <- function(b,skeleton,ylong,doselong,weight,ntreat,target,s2,j){
      Fp <- pnorm(b+qnorm(skeleton[j]))
      return(Fp*likelihood_times_prior(b,skeleton,ylong,doselong,weight,ntreat,target,s2))
      
    }
    
    } else if(model=="logit2"){ 
      
      ylong <- ylong[!is.na(ylong)]
      doselong <- doselong[!is.na(doselong)]
      weight <- weight[!is.na(weight)]
      
      ystring <- ifelse(ylong==0,yes = "N",no="T")
      ytrialr <- paste(doselong,ystring,sep="",collapse=" ") # this is the t vector for stan_crm
      
      emp_fit_rstan <- stan_crm(outcome_str=ytrialr,
                                skeleton=skeleton,
                                target=.25,
                                alpha_mean = 3,
                                alpha_sd = 10,
                                beta_mean = 0,
                                beta_sd= sqrt(1.34),
                                weights=weight,
                                model="logistic2",
                                seed=123)
                                

    }
  
  
  
  
  if(model!="logit2"){
    
    likelihood_times_prior_times_b <- function(b,skeleton,ylong,doselong,weight,ntreat,target,s2){
      return(b*likelihood_times_prior(b,skeleton,ylong,doselong,weight,ntreat,target,s2))
      
    }
    
    
    Llim <- -Inf; Ulim <- Inf;
    if(model=="prod_of_beta"){Llim <- log(1-skeleton[1]); Ulim <- Inf;}
    
    
    marginal <- integrate(likelihood_times_prior,lower = Llim,upper=Ulim,skeleton,ylong,doselong,weight,ntreat,target,s2)$value
    ndose <- length(skeleton)
    
    if(formal_bayesian){
      num_pest <- NULL
      for(j in 1:ndose){
        num_pest[j] <- integrate(likelihood_times_prior_mean,lower = Llim,upper=Ulim,skeleton,ylong,doselong,weight,ntreat,target,s2,j)$value
      }
      phat_est <- num_pest/marginal
    }
    
    
    if(!formal_bayesian){
      num_bhat <- integrate(likelihood_times_prior_times_b,lower = Llim,upper=Ulim,skeleton,ylong,doselong,weight,ntreat,target,s2)$value
      bhat_est <- num_bhat/marginal
      phat_est <- skeleton^(exp(bhat_est))
    }
    
    model_mtd <- which.min(abs(phat_est-target))

  } else{
    
    phat_est = emp_fit_rstan$prob_tox
    model_mtd = emp_fit_rstan$recommended_dose
  }
  
  
    
    
  
  
  
  
  return(list(phat_est=phat_est,model_mtd=model_mtd))
}



