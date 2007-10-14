
llHIE <- function(estimate,
                  iCoef=c(11,12,13),
                  spread=1,
                  naSpread=0.001,
                  resolution=13,
                  ...
                  ) {
   library(micEcon)
   start <- coef(estimate)
   st <- start[iCoef]
   stdds <- spread*sqrt(diag(-solve(hessian(estimate))))[iCoef]
   stdds[is.na(stdds)] <- naSpread
   var1 <- seq(st[1] - stdds[1], st[1] + stdds[1], length=resolution)
   var2 <- seq(st[2] - stdds[2], st[2] + stdds[2], length=resolution)
   var3 <- seq(st[3] - stdds[3], st[3] + stdds[3], length=resolution)
   print(rbind(var1, var2, var3))
   loglik <- array(0, dim=c(length(var1), length(var2), length(var3)))
   dimnames(loglik) <- list(var1=var1, var2=var2, var3=var3)
   data( Mroz87 )
   Mroz87$kids <- ( Mroz87$kids5 + Mroz87$kids618 > 0 )
   for(iVar1 in seq(along=var1)) {
      for(iVar2 in seq(along=var2)) {
         for(iVar3 in seq(along=var3)) {
            st <- start
            st[iCoef[1]] <- var1[iVar1]
            st[iCoef[2]] <- var2[iVar2]
            st[iCoef[3]] <- var3[iVar3]
            loglik[iVar1,iVar2,iVar3] <- selection1( lfp ~ age + I( age^2 ) + faminc + kids + educ,
                                                   wage ~ exper + I( exper^2 ) + educ + city,
                                                   data = Mroz87, start=st,
                                                   ...)
         }
      }
   }
   attr(loglik, "max") <- start[iCoef]
   attr(loglik, "grad") <- estimate$gradient[iCoef]
   class(loglik) <- c("ll3", class(loglik))
   loglik
}

compDer <- function(estimate, ...
                  ) {
   start <- coef(estimate)
   data( Mroz87 )
   Mroz87$kids <- ( Mroz87$kids5 + Mroz87$kids618 > 0 )
   selectionCD( lfp ~ age + I( age^2 ) + faminc + kids + educ,
              wage ~ exper + I( exper^2 ) + educ + city,
              data = Mroz87, start=start,
              ...)
}

plotLL3d <- function(x,
                     levels=seq(range(x, na.rm=TRUE)[1], range(x, na.rm=TRUE)[2], length=9)[-c(1,9)],
                     radius=0.01, length=2*radius) {
   library(misc3d)
   library(micEcon)
   col <- topo.colors(length(levels))
   alpha <- seq(0.2, 0.9, length=length(levels))
   cat("Levels:\n")
   print(levels)
   xVal <- as.numeric(dimnames(x)[[1]])
   yVal <- as.numeric(dimnames(x)[[2]])
   zVal <- as.numeric(dimnames(x)[[3]])
   contour3d(x,
             x=xVal, y=yVal, z=zVal,
             level=levels, color=col, alpha=alpha,
             xlab="var1", ylab="var2", zlab="var3",
             scale=TRUE
             )
   xVal <- attr(x, "max")[1]
   yVal <- attr(x, "max")[2]
   zVal <- attr(x, "max")[3]
   cat("Estimate:\n")
   print(attr(x, "max"))
   plot3d(xVal, yVal, zVal, type="s", radius=radius, col=2, alpha=1, add=T)
   grad <- attr(x, "grad")
   norm <- sqrt(sum(grad^2))
   cat("Gradient: norm =", norm, "\n")
   print(grad)
   x1 <- xVal + grad[1]/norm*length
   y1 <- yVal + grad[2]/norm*length
   z1 <- zVal + grad[3]/norm*length
   lines3d(c(xVal, x1), c(yVal, y1), c(zVal, z1), alpha=1, col=2)
   axes3d(alpha=1, color=1, labels=T)
   aspect3d(1)
}

tstHIE <- function(start=NULL, maxMethod="nr", ...) {
   cat(" ");
   if(is.null(start))
       start <- c(0,0,0,0,0,0,0,0,0,0,0,0.5,-0.5)
   library(micEcon)
   data( Mroz87 )
   Mroz87$kids <- ( Mroz87$kids5 + Mroz87$kids618 > 0 )
   mr <- selection( lfp ~ age + I( age^2 ) + faminc + kids + educ,
                   wage ~ exper + I( exper^2 ) + educ + city,
                   data = Mroz87, start=start,
                   maxMethod=maxMethod, ...)
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
   library(micEcon)
   data(Mroz87)
   Mroz87$kids <- ( Mroz87$kids5 + Mroz87$kids618 > 0 )
   mr <- selection( lfp ~ age + I( age^2 ) + faminc + kids + educ,
                   wage ~ exper + I( exper^2 ) + educ + city,
                   data = Mroz87, start = c(0,0,0,0,0,0,0,0,0,0,0,0.5,-0.5),
                   iterlim=100)
   print(summary(mr))
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
