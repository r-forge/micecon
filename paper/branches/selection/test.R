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
   curve(pnorm, -5, 5,
         ylab=expression(1 - Phi(-x^S)), xlab=expression(x^S))
   print(summary(selection(ys ~ xs, yo ~ xs)))
}

tst5 <- function(N=500, print.level=0) {
   library(mvtnorm)
   vc <- diag(3)
   vc[lower.tri(vc)] <- c(0.9, 0.5, 0.1)
   vc[upper.tri(vc)] <- vc[lower.tri(vc)]
   eps <- rmvnorm(N, rep(0, 3), vc)
   xs <- runif(N)
   ys <- xs + eps[,1] > 0
   xo1 <- runif(N)
#   xo1 <- xs
   yo1 <- xo1 + eps[,2]
   xo2 <- runif(N)
#   xo2 <- xs
   yo2 <- xo2 + eps[,3]
   a <- selection(ys~xs, list(yo1 ~ xo1, yo2 ~ xo2), print.level=print.level)
   invisible(a)
}

tstChi <- function(N=1000, rnd=0, ...) {
   library(mvtnorm)
   library(micEcon)
   set.seed(rnd)
   vc <- diag(3)
   vc[lower.tri(vc)] <- c(0.9, 0.5, 0.1)
   vc[upper.tri(vc)] <- vc[lower.tri(vc)]
   eps <- rmvnorm(N, rep(0, 3), vc)
   eps <- eps^2 - 1
   print(colMeans(eps))
   print(cor(eps))
   xs <- runif(N, -1, 1)
   ys <- xs + eps[,1] > 0
   xo1 <- runif(N)
#   xo1 <- xs
   yo1 <- xo1 + eps[,2]
   xo2 <- runif(N)
#   xo2 <- xs
   yo2 <- xo2 + eps[,3]
   a <- selection(ys~xs, list(yo1 ~ xo1, yo2 ~ xo2), ...)
   invisible(a)
}

tstChiOLS <- function(N=500, ...) {
   library(mvtnorm)
   vc <- diag(3)
   vc[lower.tri(vc)] <- c(0.9, 0.5, 0.1)
   vc[upper.tri(vc)] <- vc[lower.tri(vc)]
   eps <- rmvnorm(N, rep(0, 3), vc)
   eps <- eps^2 - 1
   print(colMeans(eps))
   print(cor(eps))
   xs <- runif(N, -1, 1)
   ys <- xs + eps[,1] > 0
   xo1 <- xs
   yo1 <- xo1 + eps[,2]
   xo2 <- runif(N)
   yo2 <- xo2 + eps[,3]
   print(summary(lm(yo1~xs, subset=ys==0)))
   print(summary(lm(yo2~xs, subset=ys==1)))
}

tstOne <- function(N=1000, print.level=0) {
   x <- runif(N, -2, 2)
   e <- rnorm(N)
   z <- x + e > 0
   y <- (x - e)*z
   a <- selection(z ~ x, y ~ x, method="2step")
   print(summary(a))
}

e1 <- function(N=200000, alpha=0.6) {
   library(mvtnorm)
   Sigma <- matrix(c(1.2,-0.5,-0.5,0.7), 2, 2)
   X <- rmvnorm(N, c(0,0), Sigma)
   x <- X[,1]
   y <- X[,2]
   sy <- sqrt(Sigma[2,2])
   ##
   cat("-- Y <", alpha, "--\n")
   my2 <- mean(y[y>alpha]^2)
   cat("MC: E[Y^2|Y >", alpha, "]=", my2, "\n")
   EY2alpha <- (sy*alpha*dnorm(alpha/sy) +
                sy^2*(1 - pnorm(alpha/sy)))/
                    (1 - pnorm(alpha/sy))
   cat("t: E[Y^2|Y >", alpha, "]=", EY2alpha, "\n")
   ##
   m <- mean(x[y > alpha]^2)
   cat("MC: E[X^2|Y >", alpha, "]=", m, "\n")
   varE <- Sigma[1,1] - Sigma[1,2]^2/sy^2
   t <- Sigma[1,2]^2/sy^4*EY2alpha + varE
   cat("t: E[X^2|Y >", alpha, "]=", t, "\n")
   ##
   cat("-- Y <", alpha, "--\n")
   my2 <- mean(y[y < alpha]^2)
   cat("MC: E[Y^2|Y <", alpha, "]=", my2, "\n")
   EY2alpha <- (-sy*alpha*dnorm(alpha/sy) +
                sy^2*pnorm(alpha/sy))/
                    pnorm(alpha/sy)
   cat("t: E[Y^2|Y <", alpha, "]=", EY2alpha, "\n")
   ##
   m <- mean(x[y < alpha]^2)
   cat("MC: E[X^2|Y <", alpha, "]=", m, "\n")
   varE <- Sigma[1,1] - Sigma[1,2]^2/sy^2
   t <- Sigma[1,2]^2/sy^4*EY2alpha + varE
   cat("t: E[X^2|Y <", alpha, "]=", t, "\n")
   ##
   cat("-- Y^2 <", alpha, "--\n")
   my2 <- mean(y[y>alpha]^2)
   cat("MC: E[Y^2|Y >", alpha, "]=", my2, "\n")
   EY2alpha <- (sy*alpha*dnorm(alpha/sy) +
                sy^2*(1 - pnorm(alpha/sy)))/
                    (1 - pnorm(alpha/sy))
   cat("t: E[Y^2|Y >", alpha, "]=", EY2alpha, "\n")
   ##
   m <- mean(x[y > alpha]^2)
   cat("MC: E[X^2|Y >", alpha, "]=", m, "\n")
   varE <- Sigma[1,1] - Sigma[1,2]^2/sy^2
   t <- Sigma[1,2]^2/sy^4*EY2alpha + varE
   cat("t: E[X^2|Y >", alpha, "]=", t, "\n")
   ##
   cat("-- Y^2 <", alpha, "--\n")
   my2 <- mean(y[y^2 < alpha]^2)
   cat("MC: E[Y^2|Y^2 <", alpha, "]=", my2, "\n")
   EY2alpha <- sy^2 - 2*sy*sqrt(alpha)*dnorm(sqrt(alpha)/sy)/(1 - 2*pnorm(-sqrt(alpha)/sy))
   cat("t: E[Y^2|Y^2 <", alpha, "]=", EY2alpha, "\n")
   ##
   m <- mean(x[y^2 < alpha]^2)
   cat("MC: E[X^2|Y^2 <", alpha, "]=", m, "\n")
   varE <- Sigma[1,1] - Sigma[1,2]^2/sy^2
   t <- Sigma[1,2]^2/sy^4*EY2alpha + varE
   cat("t: E[X^2|Y^2 <", alpha, "]=", t, "\n")
}

plotIVM <- function(N=100000) {
   EUlower <- function(alpha) {
      alpha <- sqrt(-alpha)
      EUalpha <- ss^2 - 2*ss*alpha*dnorm(alpha/ss)/(1 - 2*pnorm(-alpha/ss))
      s1s^2/ss^2*EUalpha + s1s^2 - s1s^2/ss^2
   }
   library(mvtnorm)
   vc <- diag(3)
   vc[lower.tri(vc)] <- c(0.9, 0.5, 0.1)
   vc[upper.tri(vc)] <- vc[lower.tri(vc)]
   ss <- sqrt(vc[1,1])
   s1s <- sqrt(vc[2,2])
   s12 <- vc[1,2]
   u <- rmvnorm(N, rep(0, 3), vc)
   eps <- u^2
   es <- eps[,1]
   e1 <- eps[,2]
   e2 <- eps[,3]
   curve(EUlower, -5, 0)
   for(alpha in seq(-5, 0, length=50)) {
      points(alpha, mean(e1[es < -alpha]))
   }
}
