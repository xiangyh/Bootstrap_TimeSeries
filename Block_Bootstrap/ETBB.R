##======================================
##   Extended Tapered Block Bootstrap
##======================================
# simulate data
n <- 120                               
series <- arima.sim(model=list(ar=0.2),n)                 
l <- 4                                 
B <- 1000          

stat.ETbt <- rep(NA,B) 
k <- floor(n/l)    # resample k block for each iteration
c <- 0.43          # c=0.43 has superior mean squared error performance
w <- rep(NA,l)
for(j in 1:l){
  t <- (j-0.5)/l
  if(t < c) w[j] <- t/c
  else if(t>=c & t<=(1-c)) w[j] <- 1
  else w[j] <- (1-t)/c
}
wl1 <- sum(w)
ft<-rep(0,n)
# vessel for the boostrapped values
for(irep in 1:B) {                   # the bootstrap loop
  S <- sample(0:(n-l), size=k)     # sample k beginpoint to define S_j in Shao(2011)
  for(t in 1:n){
    f <- 0
    for(j in 1:k) for(h in 1:l) f<-f+n/(k*wl1)*w[h]*ifelse(S[j]==t-h, 1,0)
    ft[t] <- f
  }
  series.ft <- ft*ts
  stat.ETbt[irep] <- mean(series.ft)
}


# > mean(stat.ETbt)
# [1] -0.05693451
# > sd(stat.ETbt)
# [1] 0.07920717
