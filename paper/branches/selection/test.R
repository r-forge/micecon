tst5 <- function(N=500, print.level=0) {
   library(mvtnorm)
   vc <- diag(3)
   vc[lower.tri(vc)] <- c(0.9, 0.5, 0.6)
   vc[upper.tri(vc)] <- vc[lower.tri(vc)]
   eps <- rmvnorm(N, rep(0, 3), vc)
   xs <- runif(N)
   ys <- xs + eps[,1] > 0
   xo1 <- runif(N)
   yo1 <- xo1 + eps[,2]
   xo2 <- runif(N)
   yo2 <- xo2 + eps[,3]
   a <- selection(ys~xs, list(yo1 ~ xo1, yo2 ~ xo2), print.level=print.level)
   invisible(a)
}

tstOne <- function(N=1000, print.level=0) {
   x <- runif(N, -2, 2)
   e <- rnorm(N)
   z <- x + e > 0
   y <- (x - e)*z
   a <- selection(z ~ x, y ~ x, method="2step")
   print(summary(a))
}
