##=============================
##   Circular Block Bootstrap
##=============================


# simulate data
n <- 120                               
series <- arima.sim(model=list(ar=0.2),n)                 
l <- 4                                 
B <- 1000          

stat.Cbt <- rep(NA,B) 
k <- floor(n/l)

series.y <- rep(series,2)
for(irep in 1:B) {                   # the bootstrap loop
  sbt <- rep(NA,k*l)                #  bootstrap replication
  for(i in 1:k) {            # random blocks
    beginpoint <- sample(1:n, size=1)     # sampling beginpoints
    sbt[(i-1)*l+1:l] <- series.y[beginpoint+0:(l-1)] # and copying blocks
  }                          # sbt has kl elements
  stat.Cbt[irep] <- stat.fun(sbt)  
}


# > mean(stat.Cbt)
# [1] -0.04355332
# > sd(stat.Cbt)
# [1] 0.1204954
