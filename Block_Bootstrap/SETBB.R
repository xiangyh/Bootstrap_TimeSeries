##===========================================
##  A Smoothing Bootstrap by Modifying ETBB
##===========================================

n <- 120                               
series <- arima.sim(model=list(ar=0.2),n)                 
l <- 4                                 
B <- 1000        

stat.SETbt <- rep(NA,B) 
k <- floor(n/l)    # resample k block for each iteration
c <- 0.43          # c=0.43 has superior mean squared error performance
bwh <- n^(-1/5)      # the bandwidth parameter by a selection method of Sheather and Jones(1991)
w <- rep(NA,l)
for(j in 1:l){
  t <- (j-0.5)/l
  if(t <= c) w[j] <- t/c
  else if(t>=c & t<=(1-c)) w[j] <- 1
  else w[j] <- (1-t)/c
}
wl1 <- sum(w)

for(irep in 1:B) {                   # the bootstrap loop
  S <- sample(0:(n-l), size=k)     # sample k beginpoint to define S_j in Shao(2011)
  for(t in 1:(n)){
    f <- 0
    for(j in 1:k) for(h in 1:l) f<-f+n/(k*wl1)*w[h]*ifelse(S[j]==t-h, 1,0)
    ft[t] <- f
  }
  z <- rnorm(n)*bwh
  series.smth <- series+z
  series.ft <- ft*series.smth
  stat.SETbt[irep] <- stat.fun(series.ft)
}


# > mean(stat.SETbt)
# [1] 0.08034387
# > sd(stat.SETbt)
# [1] 0.1196919
