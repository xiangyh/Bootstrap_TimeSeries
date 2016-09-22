##==========================================
##   Generalized Seasonal Block Bootstrap
##==========================================

# simulate data
n <- 120                               
series <- arima.sim(model=list(ar=0.2),n)                 
l <- 4                                 
B <- 1000    

d <- 6   # seasonal period
stat.GSBB.mu1 <- matrix(0, nrow=B, ncol=d) 
stat.GSBB.mu2 <- rep(NA,B) 
k <- floor(n/l)
w <- rep(floor(n/d),d)+c(rep(1,n-floor(n/d)*d),rep(0,d-n+floor(n/d)*d))    # w here is the number of occurence in the

for(irep in 1:B) {                   # the bootstrap loop
  sbt <- rep(NA,k*l)           # bootstrap replication
  for(t in (0:k)*l+1) {              # random blocks
    unif <- seq(t-d*floor((t-1)/d), t+d*floor((n-l-t)/d), by=d)
    beginpoint <- sample(unif, size=1)     # beginpoints
    sbt[t:(t+l-1)] <- series[beginpoint+0:(l-1)] # and copying blocks
  }                                 
  sbt <- sbt[1:n]
  for(i in 1:d) stat.GSBB.mu1[irep,i] <-  mean(sbt[i+0:(w[i]-1)*d])   # mu_i estimate
  stat.GSBB.mu2[irep] <- mean(stat.GSBB.mu1[irep,i])    # overall mean estimate
}


# > mean(stat.GSBB.mu1)
# [1] -0.06333954
# > mean(stat.GSBB.mu2)
# [1] 0.1219043
# > sd(stat.GSBB.mu1)
# [1] 0.3254451
# > sd(stat.GSBB.mu2)
# [1] 0.195807
