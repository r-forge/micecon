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

tstA <- function(N=500) {
   library(mvtnorm)
   eps <- rmvnorm(N, c(0,0), matrix(c(1,-0.7,-0.7,1), 2, 2))
   xs <- runif(N)
   ys <- xs + eps[,1] > 0
   yo <- (xs + eps[,2])*(ys > 0)
   a <- selection(ys ~ xs, yo ~ xs, steptol=1e-10)
   invisible(a)
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
   a <- selection(ys~xs, list(yo1 ~ xo1, yo2 ~ xo2), print.level=print.level, iterlim=30)
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

tstMroz <- function(R=10, method="ml") {
   bstat <- function(data, ind) {
      dat <- data[ind,]
      a <- selection( lfp ~ age + I( age^2 ) + faminc + kids + educ,
                     wage ~ exper + I( exper^2 ) + educ + city, data=dat)
      coef(a)
   }
   library(micEcon)
   data(Mroz87)
   Mroz87$kids <- ( Mroz87$kids5 + Mroz87$kids618 > 0 )
#   b <- boot(Mroz87, bstat, R)
   b <- selection( lfp ~ age + I( age^2 ) + faminc + kids + educ,
                  wage ~ exper + I( exper^2 ) + educ + city, data=Mroz87, method="ml")
   b
}

tstBoot <- function(R=10, N=500) {
   bstat <- function(data, ind) {
      dat <- data[ind,]
      a <- selection(ys~xs, yo ~xo, data=dat)
      coef(a)
   }
   library(micEcon)
   library(mvtnorm)
   library(boot)
   vc <- diag(2)
   vc[2,1] <- vc[1,2] <- -0.7
   eps <- rmvnorm(N, rep(0, 2), vc)
   xs <- runif(N)
   ys <- xs + eps[,1] > 0
   xo <- runif(N)
   yo <- (xo + eps[,2])*(ys > 0)
   data <- data.frame(xs,ys,xo,yo)
   b <- boot(data, bstat, R)
   b
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

plotIVM <- function(N=1000, res=50,
                    hats1=0.6166, hatr1=0.1459, hats2=1.8922, hatr2=-0.0063)
{
   EUlower <- function(alpha) {
      alpha[alpha >= 1] <- NA
      alpha <- sqrt(-alpha + 1)
      EUalpha <- ss^2 - 2*ss*alpha*dnorm(alpha/ss)/(1 - 2*pnorm(-alpha/ss))
      s1s^2/ss^2*EUalpha + s1^2 - s1s^2/ss^2 - 1
   }
   EUupper <- function(alpha) {
      alpha[alpha >= 1] <- 1
      alpha <- sqrt(-alpha + 1)
      EUalpha <- (ss*alpha*dnorm(alpha/ss) + ss^2*(1 - pnorm(alpha/ss)))/(1 - pnorm(alpha/ss))
      s2s^2/ss^2*EUalpha + s2^2 - s2s^2/ss^2 - 1
   }
   Nlower <- function(alpha) {
      alpha <- -alpha
      -Ns1s*dnorm(alpha)/pnorm(alpha)
   }
   Nupper <- function(alpha) {
      alpha <- -alpha
      Ns2s*dnorm(-alpha)/pnorm(-alpha)
   }
   library(mvtnorm)
   vc <- diag(3)
   vc[lower.tri(vc)] <- c(0.9, 0.5, 0.1)
   vc[upper.tri(vc)] <- vc[lower.tri(vc)]
   ss <- sqrt(vc[1,1])
   s1 <- sqrt(vc[2,2])
   s2 <- sqrt(vc[3,3])
   s1s <- vc[1,2]
   s2s <- vc[1,3]
   Ns1s <- hats1*hatr1
   Ns2s <- hats2*hatr2
   u <- rmvnorm(N, rep(0, 3), vc)
   eps <- u^2 - 1
   hatb1O <- coef(tmp)[c("XO1(Intercept)", "XO1xo1")]
   hatb2O <- coef(tmp)[c("XO2(Intercept)", "XO2xo2")]
   es <- eps[,1]
   e1 <- eps[,2]
   e2 <- eps[,3]
   xs <- seq(-5, 5, length=res)
   ey <- cbind(EUlower(xs), EUupper(xs),
               Nupper(cbind(1, xs)%*%hatb1O), Nlower(cbind(1,xs)%*%hatb2O),
                                        # estimated normal curves
               -s1s*dnorm(-xs)/pnorm(-xs), s2s*dnorm(xs)/pnorm(xs)
                                        # norimal curves using true parameters
               )
   mcy <- matrix(0, length(xs), 2)
   for(i in seq(length=nrow(mcy))) {
      mcy[i,1] <- mean(e1[es < -xs[i]])
      mcy[i,2] <- mean(e2[es > -xs[i]])
   }
   matplot(xs, ey, type="l", lty=c(1,1,2,2,3,3), col=1,
           xlab=expression(x^S), ylab="",
           ylim=c(-1.5,2.5))
   matpoints(xs, mcy, pch=c(1,2), cex=0.5, col=1)
   axis(4)
   legend(0, 2.4,
          legend=c(expression(
              paste("correct  ")*E*group("[", Epsilon^{O1}*group("|", Epsilon^S < -x^S, ""), "]")*
              paste("  (lower curve)"),
              phantom(paste("correct  "))*E*group("[", Epsilon^{O1}*group("|", Epsilon^S > -x^S, ""), "]")*
              paste("  (upper)"),
              E*group("[", Epsilon^{O1}*group("|", Epsilon^S < -bold(beta)^S*minute*bold(x)^S, ""), "]")*
              paste("  assumed normal (lower curve)"),
              E*group("[", Epsilon^{O1}*group("|", Epsilon^S > -bold(beta)^S*minute*bold(x)^S, ""), "]")*
              phantom(paste("  assumed normal "))*paste("  (upper)")
              )
                   ),
          bty="n", lty=c(1,1,2,2), col=1
          )
}
