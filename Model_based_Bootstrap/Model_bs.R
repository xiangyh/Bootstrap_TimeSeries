library(ggplot2)

##=========================
##  model-based bootstrap
##=========================

ar.sim<- function(n,phi){
  ts <- arima.sim(model=list(ar=phi),n)
  return(ts)
}
True.fun <- function(n,M=50000, phi=0.3){
  ar.mean <- NULL
  for(i in 1:M){
    ts <- ar.sim(n,phi)
    ar.mean <- c(ar.mean, mean(ts))
  }
  return(ar.mean)
}
modBB.fun <-  function(resid,n,B=2000){
  modBB.stat<-NULL
  for(i in 1:B){
    x.boot <- rep(0,n-1)
    resid.boot <- resid[sample(1:(n-1),size=n, replace = TRUE)]
    x.boot[1] <- resid.boot[1]/sqrt((1-phi^2))
    for(j in 2:n){
      x.boot[j] <- phi*x.boot[j-1]+resid.boot[j]
    }
    modBB.stat<-c(modBB.stat, mean(x.boot))
  }
  return(modBB.stat)
}

## the coverage ##
ci.var<-NULL
ci.mean<-NULL
for(i in 1:200){
  ts<-ar.sim(n=250,phi=0.3)
  ar.mod<-arima(ts, c(1,0,0))
  phi <- as.numeric(ar.mod$coef[1])
  resid <- ts[-1]-phi*ts[-n]
  ci.var <- c(ci.var,var(modBB.fun(resid,n)))
  ci.mean <- c(ci.mean,mean(ts))
}
rep <- 1:200
ci.right <- cbind.data.frame(mean=ci.mean,var=ci.var)
cover <- function(x) ifelse(((x[1]+1.96*sqrt(x[2])) <= 0) | ((x[1]-1.96*sqrt(x[2])) >= 0),0,1)
ci.right.cover <- apply(ci.right,1,cover)
