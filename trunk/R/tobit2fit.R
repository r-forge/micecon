tobit2fit <- function(YS, XS, YO, XO, init,
                      print.level=0, ...) {
### The model is as follows (Amemiya 1985):
### The latent variables are:
### YS* = XS'g + u1
### y2* = X'b2 + u2
### y3* = X'b3 + u3
### The observables are:
###      / 1  if  YS* > 0
### YS = \ 0  if  YS* <= 0
###      / y2*  if  YS = 0
### y2 = \ 0    if  YS = 1
###      / 0    if  YS = 0
### y3 = \ y3*  if  YS = 1
### 
###  YS      binary or logical vector, 0 (FALSE) corresponds to observable y2*,
###          1 (TRUE) to observable y3*
###  y2, y3  numeric vector, outcomes
###  X       explanatory variables for outcomes
###  XS              -"-                selection, should include exclusion restriction
###  ...     additional parameters for maxLik
### Result:
### Object of class 'tobit5', derived from 'maxLik'.
### Includes all the components of maxLik and additionally
### twoStep   Results for Heckman two-step estimator, used for initial values
### 
   loglik <- function( beta) {
      g <- beta[igamma]
      b <- beta[ibeta]
      sigma <- beta[isigma]
      if(sigma < 0) return(NA)
      rho <- beta[irho]
      if( ( rho < -1) || ( rho > 1)) return(NA)
      Z1.g <- Z1 %*% g
      Z2.g <- Z2 %*% g
      X2.b <- X2 %*% b
      u2 <- Y21 - X2.b
      r <- sqrt( 1 - rho^2)
      B <- (Z2.g + rho/sigma*u2)/r
      l1 <- sum(pnorm(-Z1.g, log.p=TRUE))
                                        # loglik related with unobserved cases ...
      l2 <- -N2/2*log(2*pi) - N2*log(sigma) +
          sum(pnorm(B, log.p=TRUE) - 0.5*(u2/sigma)^2)
                                        # ... and observed cases
      loglik <- l1 + l2
   }
    gradlik <- function(beta) {
       g <- beta[igamma]
       b <- beta[ibeta]
       sigma <- beta[isigma]
       if(sigma < 0) return(NA)
       rho <- beta[irho]
       if( ( rho < -1) || ( rho > 1)) return(NA)
       Z1.g <- Z1 %*% g
       Z2.g <- Z2 %*% g
       X2.b <- X2 %*% b
       u2 <- Y21 - X2.b
       r <- sqrt( 1 - rho^2)
       B <- (Z2.g + rho/sigma*u2)/r
       lambdaB <- dnorm(B)/pnorm(B)
       gradient <- numeric(NParam)
       gradient[igamma] <- t(Z1) %*% (-dnorm(-Z1.g)/pnorm(-Z1.g)) +
           (t(Z2) %*% lambdaB)/r
       gradient[ibeta] <- t(X2) %*% (u2/sigma^2 - lambdaB*rho/sigma/r)
       gradient[isigma] <- sum(u2^2/sigma^3 - lambdaB*rho*u2/sigma^2/r) - N2/sigma
       gradient[irho] <- sum(lambdaB*(u2/sigma + rho*Z2.g))/r^3
       gradient
    }
    hesslik <- function(beta) {
       g <- beta[igamma]
       b <- beta[ibeta]
       sigma <- beta[isigma]
       if(sigma < 0) return(NA)
       rho <- beta[irho]
       if( ( rho < -1) || ( rho > 1)) return(NA)
       Z1.g <- as.vector(Z1 %*% g)
       Z2.g <- as.vector(Z2 %*% g)
       X2.b <- as.vector(X2 %*% b)
       u2 <- Y21 - X2.b
       r <- sqrt( 1 - rho^2)
       B <- (Z2.g + rho/sigma*u2)/r
       lambdaB <- dnorm(B)/pnorm(B)
       fZ1.g <- dnorm(-Z1.g)
       FZ1.g <- pnorm(-Z1.g)
                                        # the previous code
                                        #       C <- ifelse(B > -25,
                                        #                   -(pnorm(B)*dnorm(B)*B + dnorm(B)^2)/pnorm(B)^2,
                                        #                   -1)
       C <- ifelse(B > -500,
                   -exp(dnorm(B, log = TRUE) - pnorm(B, log.p = TRUE))*B -
                   exp(2 * (dnorm(B, log = TRUE) - pnorm(B, log.p = TRUE))),
                   -1)
                                        # recommended by Dimitrios Rizopoulos, KULeuven
                                        # This is a hack in order to avoid numerical problems.  How to do
                                        # it better?  How to prove the limit value?
       hess <- matrix(0, NParam, NParam)
       a <- (-fZ1.g*FZ1.g*Z1.g + fZ1.g^2)/FZ1.g^2
       hess[igamma,igamma] <- -t(Z1) %*% (Z1*a) + t(Z2) %*% (Z2*C)/r^2
       hess[igamma,ibeta] <- -t(Z2) %*% (X2*C)*rho/r^2/sigma
       hess[ibeta,igamma] <- t(hess[igamma,ibeta])
       hess[igamma,isigma] <- -rho/sigma^2/r^2*t(Z2) %*% (C*u2)
       hess[isigma,igamma] <- t(hess[igamma,isigma])
       hess[igamma,irho] <- t(Z2) %*%
           (C*(u2/sigma + rho*Z2.g)/r^4 + lambdaB*rho/r^3)
       hess[irho,igamma] <- t(hess[igamma,irho])
       hess[ibeta,ibeta] <- t(X2) %*%
           (X2 * ((rho/r)^2*C - 1))/sigma^2
       hess[ibeta,isigma] <- t(X2) %*%
           (C*rho^2/sigma^3*u2/r^2 +
            rho/sigma^2*lambdaB/r - 2*u2/sigma^3)
       hess[isigma,ibeta] <- t(hess[ibeta,isigma])
       hess[ibeta,irho] <- t(X2) %*%
           (-C*(u2/sigma + rho*Z2.g)/r^4*rho -
            lambdaB/r^3)/sigma
       hess[irho,ibeta] <- t(hess[ibeta,irho])
       hess[isigma,isigma] <- sum(
                                  -3*u2*u2/sigma^4
                                  +2*lambdaB* u2/r *rho/sigma^3
                                  +rho^2/sigma^4 *u2*u2/r^2 *C) +
                                      N2/sigma^2
       hess[isigma,irho] <- hess[irho,isigma] <-
           -sum((C*rho*(u2/sigma + rho*Z2.g)/r + lambdaB)*
                u2/sigma^2)/r^3
       hess[irho,irho] <-
           sum(C*((u2/sigma + rho*Z2.g)/r^3)^2 +
               lambdaB*(Z2.g*(1 + 2*rho^2) + 3*rho*u2/sigma) / r^5 )
       return( hess )
    }
    ## ---------------
    NXS <- ncol( XS)
    if(is.null(colnames(XS)))
        colnames(XS) <- rep("XS", NXS)
    NXO <- ncol( XO)
    if(is.null(colnames(XO)))
        colnames(XO) <- rep("XO", NXO)
    Nparam <- NXS + NXO + 2
                                        # Total # of parameters
    Nobs <- length( YS)
    NO <- length( YS[YS > 0])
    ## indices in for the parameter vector
    igamma <- 1:NXS
    ibeta <- seq(tail(igamma, 1)+1, length=NXO)
    isigma <- tail(ibeta, 1) + 1
    irho <- tail(isigma, 1) + 1
    ## split data by selection
    YO <- YO[YS==0]
    XO <- XO[YS==0,,drop=FALSE]
    XS1 <- XS[YS==0,,drop=FALSE]
    XS2 <- XS[YS==1,,drop=FALSE]
    if(print.level > 0) {
        cat( "Yo observed:", NO, "times; not observed:", NObs - NO, "times\n")
    }
    if( print.level > 0) {
        cat( "Initial values:\n")
        cat("beta1\n")
        print(init[ibeta])
        cat("sigma1\n")
        print(init[isigma])
        cat("rho1\n")
        print(init[irho])
    }
   estimation <- maxLik(loglik, grad=gradlik, hess=hesslik, theta=init,
                     print.level=print.level, ...)
    result <- c(estimation,
                twoStep=probit
                                        # should be something like heckit5
                )
    class(result) <- c("tobit2", class(estimation))
    result
}
