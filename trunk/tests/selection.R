library(micEcon)
library(mvtnorm)
N <- 1500
NNA <- 5
vc <- diag(3)
vc[lower.tri(vc)] <- c(0.9, 0.5, 0.6)
vc[upper.tri(vc)] <- vc[lower.tri(vc)]
set.seed(1)
## ------- Tobit-5 example ---------
eps <- rmvnorm(N, rep(0, 3), vc)
xs <- runif(N)
ys <- xs + eps[,1] > 0
xo1 <- runif(N)
yo1 <- xo1 + eps[,2]
xo2 <- runif(N)
yo2 <- xo2 + eps[,3]
## Put some NA-s into the data
ys[sample(N, NNA)] <- NA
xs[sample(N, NNA)] <- NA
xo1[sample(N, NNA)] <- NA
xo2[sample(N, NNA)] <- NA
yo1[sample(N, NNA)] <- NA
yo2[sample(N, NNA)] <- NA
a <- selection(ys~xs, list(yo1 ~ xo1, yo2 ~ xo2))
print(summary(a))
## ------- Tobit-2 exmple -----------
vc <- diag(2)
vc[2,1] <- vc[1,2] <- -0.7
eps <- rmvnorm(N, rep(0, 2), vc)
xs <- runif(N)
ys <- xs + eps[,1] > 0
xo <- runif(N)
yo <- (xo + eps[,2])*(ys > 0)
xs[sample(N, NNA)] <- NA
ys[sample(N, NNA)] <- NA
xo[sample(N, NNA)] <- NA
yo[sample(N, NNA)] <- NA
a <- selection(ys~xs, yo ~xo)
print(summary(a))

