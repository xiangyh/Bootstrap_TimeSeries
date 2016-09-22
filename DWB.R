##======================================
##     Dependent Wild Bootstrap
##======================================

library(MASS)
# simulate data
n <- 120                               
series <- arima.sim(model=list(ar=0.2),n)                 
l <- 4                                 
B <- 1000     

stat.DWbt <- rep(NA,B) 
covmat.Bt <- matrix(0,nrow=n,ncol=n)       # covariance matrix of Wt by Bartlett kernel 
for(row in 1:n){
  for(col in 1:n){
    covmat.Bt[row, col] <- (1-abs(row-col)/l)*ifelse(abs(row-col)/l<=1, 1,0)
  }
}

c <- 0.43
wtrap <- function(x){                       # define w(x)
  x/c*ifelse((x<c)&(x>=0),1,0)+ifelse((x>=c)&(x<=(1-c)),1,0)+(1-x)/c*ifelse((x>(1-c))&(x<=1),1,0)
}
wwt <- function(t){                         # define w.w(t)
  wwt.fun <- function(x) wtrap(x)*wtrap(x+t)
  integrate(wwt.fun,lower=-1,upper=1)$value            
}

covmat.wt <- matrix(0,nrow=n,ncol=n)       # covariance matrix of Wt by self convolution of w(t)
wwt0 <- wwt(0)
for(row in 1:n){
  for(col in 1:n){
    rc <- abs(col-row)/l
    covmat.wt[row,col] <- wwt(rc)/wwt0
  }
}

stat.DWbt <- rep(NA, B)
for(irep in 1:B) {                   # the bootstrap loop
  S <- sample(0:(n-l), size=k)     # sample k beginpoint to define S_j in Shao(2011)
  axil.w <- mvrnorm(1, mu=rep(0,n), Sigma=covmat.wt)
  series.bt <- mean(series) + (series-mean(series))*axil.w  # DWB pseudo observations for each B
  stat.DWbt[irep] <- mean(series.bt)
}


# > mean(stat.DWbt)
# [1] 0.03121384
# > sd(stat.DWbt)
# [1] 0.1115432