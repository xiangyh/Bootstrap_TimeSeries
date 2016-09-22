
##=============================
##   Moving Block Bootstrap
##=============================

# simulate data
n <- 120                               
series <- arima.sim(model=list(ar=0.2),n)                 
l <- 4                                 
B <- 1000                            

stat.Mbt <- rep(NA,B) 
k <- floor(n/l) # size

# the bootstrap loop
for(irep in 1:B) {                   
  sbt <- rep(NA,k*l)     
  
  # random blocks
  for(i in 1:k) {           
    beginpoint <- sample(1:(n-l+1), size=1)     # by randomly sampling beginpoints
    sbt[(i-1)*l+1:l] <- ts[beginpoint+0:(l-1)] # and copying blocks
  }                                 # sbt has kl <= l elements.
  stat.Mbt[irep] <- mean(sbt)  # the autocorrlation estimate
}

# > mean(stat.Mbt)
# [1] -0.09734198

# > sd(stat.Mbt)
# [1] 0.1264649

# > hist(stat.Mbt)
