
# feasibility function

HFD_check <- function(X,feas_targ,pu_feas,prior_feas,aplus){
  
  # X is the highest feasible dose count. 
  # feas_targ is the minimum required feasibility probability
  # pu_feas is the upper probabiltiy cuttoff for feasibility
  # prior feas is the prior values for each dose being feasible - q0 
  
  ndose <- length(X)-1
  doses <- 1:ndose
  #q0 <- c(1,.975,.95,.90,.75,.50) # note no 0 at the end here 
  pi0 <- c(-1*diff(prior_feas),prior_feas[length(prior_feas)])
  
  aplus <- 2
  a <- pi0*aplus
  alf <- bet <- ptheta <- NULL
  for(j in 1:ndose){
    alf[j] <- sum(a[(j+1):(ndose+1)])+sum(X[(j+1):(ndose+1)])
    bet[j] <- sum(a[1:j]) + sum(X[1:j])
    ptheta[j] <- pbeta(q=feas_targ,shape1 = alf[j],shape2 = bet[j])
  }
  
  HFD <- sum(ptheta<pu_feas)
  
  return(HFD)
  
}

  
  
  