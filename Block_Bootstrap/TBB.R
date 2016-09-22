
##=============================
##   Tapered Block Bootstrap
##=============================  

# simulate data
n <- 120                               
series <- arima.sim(model=list(ar=0.2),n)                 
l <- 4                                 
B <- 1000                            

stat.Tbt <- rep(NA,B) 
k <- floor(n/l)    # resample k block for each iteration
c <- 0.43          # c=0.43 has superior mean squared error performance
w <- rep(NA,l)
for(j in 1:l){
  t <- (j-0.5)/l
  if(t <= c) w[j] <- t/c
  else if(t>=c & t<=(1-c)) w[j] <- 1
  else w[j] <- (1-t)/c
}
f <- w*l^0.5/(sum(w^2)^0.5)

for(irep in 1:B) {                   # the bootstrap loop
  series.bt <- rep(NA,k*l)                # bootstrap replication
  for(i in 1:k) {            # random blocks
    beginpoint <- sample(1:(n-l+1), size=1)     # sampling beginpoints
    series.bt[(i-1)*l+1:l] <- f*series[beginpoint+0:(l-1)] # and copying blocks
  }
  stat.Tbt[irep] <- stat.fun(series.bt) 
}

# > mean(stat.Tbt)
# [1] 0.03788246
# > sd(stat.Tbt)
# [1] 0.09401496
