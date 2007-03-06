tst2 <- function(N=1000, rho=-0.7, print.level=0) {
   library(mvtnorm)
   vc <- diag(2)
   vc[2,1] <- vc[1,2] <- -0.7
   eps <- rmvnorm(N, rep(0, 2), vc)
   xs <- runif(N)
   ys <- xs + eps[,1] > 0
   xo <- runif(N)
   yo <- (xo + eps[,2])*(ys > 0)
   a <- selection(ys~xs, yo ~xo, print.level=print.level)
   invisible(a)
}

tstA <- function() {
   library(mvtnorm)
   eps <- rmvnorm(500, c(0,0), matrix(c(1,-0.7,-0.7,1), 2, 2))
   xs <- runif(500)
   ys <- xs + eps[,1] > 0
   yo <- (xs + eps[,2])*(ys > 0)
   print(summary(selection(ys ~ xs, yo ~ xs)))
}

tstB <- function() {
   library(mvtnorm)
   eps <- rmvnorm(500, c(0,0), matrix(c(1,-0.7,-0.7,1), 2, 2))
   xs <- runif(500, -5, 5)
   ys <- xs + eps[,1] > 0
   yo <- (xs + eps[,2])*(ys > 0)
   pdf(file="pnorm.pdf", width=6, height=4)
   curve(pnorm, -5, 5,
         ylab=expression(1 - Phi(-x^S)), xlab=expression(x^S))
   dev.off()
   print(summary(selection(ys ~ xs, yo ~ xs)))
}

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
