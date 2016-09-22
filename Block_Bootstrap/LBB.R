##===========================
##  Local Block Bootstrap
##===========================

# simulate data
n <- 120                               
series <- arima.sim(model=list(ar=0.2),n)                 
l <- 4                                 
B <- 1000    

stat.Lbt <- matrix(0,nrow=B, ncol=n) 
k <- floor(n/l)
prob <- n^(-0.3)
h <- n^(-0.3)
edge <- floor(n*prob)
# vessel for the boostrapped values
for(irep in 1:B) {                   # the bootstrap loop
  series.bt <- rep(NA,k*l)                # local vector for a bootstrap replication
  for(i in 1:k) {            # fill the vector with random blocks
    beginpoint <- sample(-edge:edge, size=1)     # by randomly sampling endpoints
    # if j+(i-1)b+Ii is outside the range of interegrs 1 to n
    for(j in 1:l){
      index <- beginpoint + (i-1)*l + j
      if(index<1|index>n) series.bt[(i-1)*l+j] <- series[-beginpoint + (i-1)*l + j]      
      else series.bt[(i-1)*l+j] <- series[index]
    }
  }
  for(t in 1:n){
    x <- (series[t]-series)/h
    x <- ifelse(x<--0.5|x>0.5,0.5,x)
    kx <- 1/6*(x+0.5)*(x+0.5)
    stat.Lbt[irep,t] <- sum(kx*series.bt)/(n*h)
  }
}

# > mean(stat.Lbt)
# [1] -0.04131238
# > sd(stat.Lbt)
# [1] 0.07297704

