# dissertation curve generation functions

# conaway curves and bagley curves


conawayclass <- function(ndoses,maxtox,mintox,numcurves){
  x=seq(0.001,0.999,0.001)
  f=function(x,alpha,beta,kappa){
    beta+((alpha-beta)/(1+exp(-kappa*(log(x/(1-x))))))
  }
  lb=ub=xstar=fstar=rep(0,ndoses)
  for(i in 1:ndoses){
    lb[i]=ceiling((i-1)*(999/ndoses))
    ub[i]=floor(i*(999/ndoses))
  }
  lb[1]=1
  curves=matrix(,nrow=numcurves,ncol=ndoses)
  for(k in 1:numcurves){
    beta=runif(1,0,mintox)
    alpha=runif(1,beta,maxtox)
    kappa=runif(1,1,2)
    for(j in 1:ndoses){
      xstar[j]=x[sample(lb[j]:ub[j],1)]
      fstar[j]=f(xstar[j],alpha,beta,kappa)
    }
    curves[k,]=fstar
  }
  curves
} # generates toxicity curves
curve_gen_2 <- function(ndose,ncurve,b){
  U <- runif(n=ncurve,min=0,max=b)
  
  curves <- matrix(NA,ncol = ndose,nrow=ncurve)
  
  for(i in 1:ncurve){
    for(j in 1:ndose){
      if(j==1){
        curves[i,j] <- runif(n=1,min=1-U[i],max=1)
      } else{
        curves[i,j] <- runif(n=1,min=curves[i,j-1]-U[i],max=curves[i,j-1])
      }
    }
    curves[i,][curves[i,]<0] <- .001
    
  }
  
  
  return(curves)
  
}