piecewise.generate <-
  function(N, het="discrete",
           NX=3,
           cens=TRUE, cens.time=c(5,20)) {
### N   - number of obs. to generate
### het = c("discrete", "gamma", "none")
###       type of heterogeneity
### NX  - # of explanatory variables
### cens - whether to censore observations
### cens.time - when to censore.  It must be a vector of length 2, censoring
###        times are given as runif(N, cens.time[1], cens.time[2])
###
### Return:
### A list with the following components:
### X    - NxNX matrix of explanatory variables
### v    - N-vector of heterogeneity component
### t    - vector of durations
### cens - logical vector whether observation is censored
    gamma <- c(1,-1,0.5,-0.5,1.5,-1.5)
    X <- matrix(runif(N*NX), N, NX)
    v <- switch(het,
              discrete = -1.5*rbinom(N,1,0.3),
              gamma = log(rgamma(N, 4, scale=1/4)),
              none = 0)
  mu <- exp(v + X %*% gamma[1:NX])
  t <- rexp(N, 0.2*mu)
  t[t>2] <- 2 + rexp(N, 0.1*mu)[t>2]
  t[t>5] <- 5 + rexp(N, 0.05*mu)[t>5]
  censored <- rep(0, N)
  if(cens) {
    if(length(cens.time) == 1)
      tc <- rep(cens.time, N)
    else {
      tc <- runif(N,cens.time[1], cens.time[2])
    }
    censored[t>tc] <- 1
    t[t>tc] <- tc[t>tc]
  }
  invisible(list(X=X,
       v=v,
       t=t,
       cens=censored))
}

