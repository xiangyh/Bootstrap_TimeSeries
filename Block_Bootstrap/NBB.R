##===================================
##   Nonoverlapping Block Bootstrap
##===================================

# simulate data
n <- 120                               
series <- arima.sim(model=list(ar=0.2),n)                 
l <- 4                                 
B <- 1000        
stat.fun<-function(x){
  mean(x)
}

stat.Nbt <- rep(NA,B)  
b <- floor(n/l)
# vessel for the boostrapped values
for(irep in 1:B) {                   # the bootstrap loop
  sbt <- rep(NA,b*l)                #  bootstrap replication
  for(i in 1:b) {            # random blocks
    subpoint <- (0:(b-1))*l+1
    beginpoint <- sample(subpoint, size=1)     # sampling beginpoints
    sbt[(i-1)*l+1:l] <- series[beginpoint+0:(l-1)] # and copying blocks
  }
  stat.Nbt[irep] <- stat.fun(sbt)  # the target_stat. estimate
}


# > mean(stat.Nbt)
# [1] -0.06902236
# > sd(stat.Nbt)
# [1] 0.07680826

