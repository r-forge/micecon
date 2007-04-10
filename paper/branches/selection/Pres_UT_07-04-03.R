ex1 <- function(N=500) {
   ## tobit2, correct specification, exclusion restriction
   library(micEcon)
   library(mvtnorm)
   cor <- -0.9
   eps <- rmvnorm(N, c(0,0), matrix(c(1, cor, cor, 1), 2, 2))
   xs <- runif(N)
   ys <- xs + eps[,1] > 0
   xo <- runif(N)
   yoX <- xo + eps[,2]
   yo <- yoX*(ys > 0)
   print(summary(m <- selection(ys~xs, yo ~xo)))
##    print(summary(lm(yo ~ xo, subset=ys == 1)))
   col <- c("lightgray", "black")
   plot(xo, yoX, col=col[1 + ys], cex=0.5)
   abline(a=0, b=1, lwd=2)
   abline(a=coef(m)[3], b=coef(m)[4], col=3, lwd=2)
   cf <- coef(lm(yo ~ xo, subset=ys==1))
   abline(a=cf[1], b=cf[2], col=2, lwd=2)
}

ex2 <- function(N=1500) {
   ## tobit2, correct specification, no exclusion restriction
   library(micEcon)
   library(mvtnorm)
   cor <- -0.7
   eps <- rmvnorm(N, c(0,0), matrix(c(1, cor, cor, 1), 2, 2))
   xs <- runif(N)
   ys <- xs + eps[,1] > 0
   xo <- xs
   yoX <- xo + eps[,2]
   yo <- yoX*(ys > 0)
   print(summary(m <- selection(ys~xs, yo ~xo)))
##    print(summary(lm(yo ~ xo, subset=ys == 1)))
   col <- c("lightgray", "black")
   plot(xo, yoX, col=col[1 + ys], cex=0.5)
   abline(a=0, b=1, lwd=2)
   abline(a=coef(m)[3], b=coef(m)[4], col=3, lwd=2)
   cf <- coef(lm(yo ~ xo, subset=ys==1))
   abline(a=cf[1], b=cf[2], col=2, lwd=2)
}   

ex2a <- function(N=500) {
   ## tobit2, correct specification, no exclusion restriction
   library(micEcon)
   library(mvtnorm)
   cor <- -0.9
   eps <- rmvnorm(N, c(0,0), matrix(c(1, cor, cor, 1), 2, 2))
   xs <- runif(N, -5, 5)
   ys <- xs + eps[,1] > 0
   xo <- xs
   yoX <- xo + eps[,2]
   yo <- yoX*(ys > 0)
   print(summary(m <- selection(ys~xs, yo ~xo)))
##    print(summary(lm(yo ~ xo, subset=ys == 1)))
   col <- c("lightgray", "black")
   plot(xo, yoX, col=col[1 + ys], cex=0.5)
   abline(a=0, b=1, lwd=2)
   abline(a=coef(m)[3], b=coef(m)[4], col=3, lwd=2)
   cf <- coef(lm(yo ~ xo, subset=ys==1))
   abline(a=cf[1], b=cf[2], col=2, lwd=2)
}   

ex3 <- function(N=1500) {
   ## tobit2, misspecification, no exclusion restriction
   library(micEcon)
   library(mvtnorm)
#   set.seed(16)
   cor <- 0.9
   eps <- rmvnorm(N, c(0,0), matrix(c(1, cor, cor, 1), 2, 2))
   eps <- eps^2 - 1
                                        # chi2(1) distribution with mean zero
   xs <- runif(N, -1, 1)
   ys <- xs + eps[,1] > 0
   xo <- runif(N)
   yoX <- xo + eps[,2]
   yo <- yoX*(ys > 0)
   print(summary(m <- selection(ys~xs, yo ~xo)))
##    print(summary(lm(yo ~ xo, subset=ys == 1)))
   col <- c("lightgray", "black")
   plot(xo, yoX, col=col[1 + ys], cex=0.5)
   abline(a=0, b=1, lwd=2)
   abline(a=coef(m)[3], b=coef(m)[4], col=3, lwd=2)
   cf <- coef(lm(yo ~ xo, subset=ys==1))
   abline(a=cf[1], b=cf[2], col=2, lwd=2)
}   

ex4 <- function(N=1500) {
   ## tobit2, correct specification, no exclusion restriction
   library(micEcon)
   library(mvtnorm)
#   set.seed(16)
   cor <- -0.9
   eps <- rmvnorm(N, c(0,0), matrix(c(1, cor, cor, 1), 2, 2))
   eps <- eps^2 - 1
                                        # chi2(1) distribution with mean zero
   xs <- runif(N, -1, 1)
   ys <- xs + eps[,1] > 0
   xo <- xs
   yoX <- xo + eps[,2]
   yo <- yoX*(ys > 0)
   print(summary(m <- selection(ys~xs, yo ~xo, iterlim=30)))
##    print(summary(lm(yo ~ xo, subset=ys == 1)))
   col <- c("lightgray", "black")
   plot(xo, yoX, col=col[1 + ys], cex=0.5)
   abline(a=0, b=1, lwd=2)
   abline(a=coef(m)[3], b=coef(m)[4], col=3, lwd=2)
   cf <- coef(lm(yo ~ xo, subset=ys==1))
   abline(a=cf[1], b=cf[2], col=2, lty=2, lwd=2)
}   
