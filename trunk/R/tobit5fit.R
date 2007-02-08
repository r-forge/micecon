tobit5fit <- function(YS, XS, YO1, XO1, YO2, XO2, start,
                      print.level=0, ...) {
### The model is as follows (Amemiya 1985):
### The latent variables are:
### YS* = XS'g + u1
### y2* = X'b1 + u2
### y3* = X'b2 + u3
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
        b1 <- beta[ibeta1]
        sigma2 <- beta[isigma1]
        rho2 <- beta[irho1]
        if( ( rho2 < -1) || ( rho2 > 1)) return(NA)
        b2 <- beta[ibeta2]
        sigma3 <- beta[isigma2]
        rho3 <- beta[irho2]
        if((rho3 < -1) || (rho3 > 1)) return(NA)
                                        # check the range
        Z2g <- XS1%*%g
        Z3g <- XS2%*%g
        XO1b1 <- XO1%*%b1
        XO2b2 <- XO2%*%b2
        u2 <- YO1 - XO1b1
        u3 <- YO2 - XO2b2
        sqrt1r22 <- sqrt( 1 - rho2^2)
        sqrt1r32 <- sqrt( 1 - rho3^2)
        B1 <- -(Z2g + rho2/sigma2*u2)/sqrt1r22
        B2 <- (Z3g + rho3/sigma3*u3)/sqrt1r32
    l2 <- sum(
              -log( sigma2) - 0.5*( u2/sigma2)^2 +
              log( pnorm( B1)))
    l3 <- sum(
              -log( sigma3) - 0.5*( u3/sigma3)^2 +
              log( pnorm( B2)))
    loglik <- l2 + l3 - Nobs/2*log( 2*pi)
}
    gradlik <- function(beta) {
       ## gradient is 1xNParam matrix
        ## components of the gradient are ordered as: g b_2, s_2, r_2, b_3,
        ## s_3, r_3
        g <- beta[igamma]
        b1 <- beta[ibeta1]
        sigma2 <- beta[isigma1]
        rho2 <- beta[irho1]
        b2 <- beta[ibeta2]
        sigma3 <- beta[isigma2]
        rho3 <- beta[irho2]
        Z2g <- XS1%*%g
        Z3g <- XS2%*%g
        XO1b1 <- XO1%*%b1
        XO2b2 <- XO2%*%b2
        u2 <- YO1 - XO1b1
        u3 <- YO2 - XO2b2
        sqrt1r22 <- sqrt( 1 - rho2^2)
        sqrt1r32 <- sqrt( 1 - rho3^2)
        B1 <- -(Z2g + rho2/sigma2*u2)/sqrt1r22
        B2 <- (Z3g + rho3/sigma3*u3)/sqrt1r32
        fB1 <- dnorm( B1)
        FB1 <- pnorm( B1)
        fB2 <- dnorm( B2)
        FB2 <- pnorm( B2)
        lambda1 <- fB1/FB1
        lambda2 <- fB2/FB2
                                        # now the calculation
        l.g <- -t(lambda1) %*% XS1/sqrt1r22 + t(lambda2) %*% XS2/sqrt1r32
        l.b1 <-  t(lambda1*rho2/sigma2/sqrt1r22 + u2/sigma2^2) %*% XO1
        l.b2 <- t(-lambda2*rho3/sigma3/sqrt1r32 + u3/sigma3^2) %*% XO2
        l.s2 <- sum( -1/sigma2 + u2^2/sigma2^3
                    +lambda1*rho2/sigma2^2*u2/sqrt1r22)
        l.s3 <- sum( -1/sigma3 + u3^2/sigma3^3
                    -lambda2*rho3/sigma3^2*u3/sqrt1r32)
        l.r2 <- -sum( lambda1*( u2/sigma2 + rho2*Z2g)/sqrt1r22^3)
        l.r3 <- sum( lambda2*( u3/sigma3 + rho3*Z3g)/sqrt1r32^3)
        gradient <- cbind( l.g, l.b1, l.s2, l.r2, l.b2, l.s3, l.r3)
        gradient
    }
    hesslik <- function(beta) {
        g <- beta[igamma]
        b1 <- beta[ibeta1]
        sigma2 <- beta[isigma1]
        rho2 <- beta[irho1]
        b2 <- beta[ibeta2]
        sigma3 <- beta[isigma2]
        rho3 <- beta[irho2]
        Z2g <- XS1%*%g
        Z3g <- XS2%*%g
        XO1b1 <- XO1%*%b1
        XO2b2 <- XO2%*%b2
        u2 <- YO1 - XO1b1
        u3 <- YO2 - XO2b2
        sqrt1r22 <- sqrt( 1 - rho2^2)
        sqrt1r32 <- sqrt( 1 - rho3^2)
        B1 <- -(Z2g + rho2/sigma2*u2)/sqrt1r22
        B2 <- (Z3g + rho3/sigma3*u3)/sqrt1r32
        fB1 <- dnorm( B1)
        FB1 <- pnorm( B1)
        fB2 <- dnorm( B2)
        FB2 <- pnorm( B2)
        lambda1 <- fB1/FB1
        lambda2 <- fB2/FB2
        CB1 <- as.vector( -(FB1*fB1*B1 + fB1*fB1)/FB1/FB1)
        CB2 <- as.vector( -(FB2*fB2*B2 + fB2*fB2)/FB2/FB2)
                                        # now the calculation
        l.gg <- t( XS1) %*% ( XS1 * CB1)/sqrt1r22^2 +
            t( XS2) %*% ( XS2 * CB2)/sqrt1r32^2
        l.gb1 <- -t( XS1) %*%
            ( XO1 * CB1)*rho2/sqrt1r22^2/sigma2
        l.gs2 <- -rho2/sigma2^2/sqrt1r22^2*
            t( XS1) %*% ( CB1*u2)
        l.gr2 <- t( XS1) %*%
            ( CB1*( u2/sigma2 + rho2*Z2g)/sqrt1r22^4 -
             lambda1*rho2/sqrt1r22^3)
        l.gb2 <- -t( XS2) %*%
            ( XO2 * CB2)*rho3/sqrt1r32^2/sigma3
        l.gs3 <- -rho3/sigma3^2/sqrt1r32^2*
            t( XS2) %*% ( CB2*u3)
        l.gr3 <- t( XS2) %*%
            ( CB2*( u3/sigma3 + rho3*Z3g)/sqrt1r32^4 +
             lambda2*rho3/sqrt1r32^3)
        l.b1b1 <- t( XO1) %*%
            (XO1 * ( (rho2/sqrt1r22)^2 * CB1 - 1))/sigma2^2
        l.b1s2 <- t( XO1) %*%
            ( CB1*rho2^2/sigma2^3*u2/sqrt1r22^2 -
             rho2/sigma2^2*lambda1/sqrt1r22 -
             2*u2/sigma2^3)
        l.b1r2 <- t( XO1) %*%
            ( -CB1*( u2/sigma2 + rho2*Z2g)/sqrt1r22^4*rho2 +
             lambda1/sqrt1r22^3)/sigma2
        ## l.b1x3 is zero
        l.s2s2 <- sum(
                      1/sigma2^2
                      -3*u2*u2/sigma2^4
                      + u2*u2/sigma2^4 *rho2^2/sqrt1r22^2 *CB1
                      -2*lambda1* u2/sqrt1r22 *rho2/sigma2^3)
        l.s2r2 <- sum(
                      ( -CB1*rho2*(u2/sigma2 + rho2*Z2g)/sqrt1r22 +
                       lambda1)
                      *u2/sigma2^2)/sqrt1r22^3
        ## l.s2x3 is zero
        l.r2r2 <- sum(
                      CB1*( ( u2/sigma2 + rho2*Z2g)/sqrt1r22^3)^2
                      -lambda1*( Z2g*( 1 + 2*rho2^2) + 3*rho2*u2/sigma2) /
                      sqrt1r22^5
                      )
        ## l.r2x3 is zero
        l.b2b2 <- t( XO2) %*%
            (XO2 * ( (rho3/sqrt1r32)^2 * CB2 - 1))/sigma3^2
        l.b2s3 <- t( XO2) %*%
            ( CB2*rho3^2/sigma3^3*u3/sqrt1r32^2 +
             rho3/sigma3^2*lambda2/sqrt1r32 - 2*u3/sigma3^3)
        l.b2r3 <- t( XO2) %*%
            ( -CB2*( u3/sigma3 + rho3*Z3g)/sqrt1r32^4*rho3 -
             lambda2/sqrt1r32^3)/sigma3
        l.s3s3 <- sum(
                      1/sigma3^2
                      -3*u3*u3/sigma3^4
                      +2*lambda2* u3/sqrt1r32 *rho3/sigma3^3
                      +rho3^2/sigma3^4 *u3*u3/sqrt1r32^2 *CB2)
        l.s3r3 <- -sum(
                       ( CB2*rho3*(u3/sigma3 + rho3*Z3g)/sqrt1r32 +
                        lambda2)
                       *u3/sigma3^2)/sqrt1r32^3
        l.r3r3 <- sum(
                      CB2*( ( u3/sigma3 + rho3*Z3g)/sqrt1r32^3)^2
                      + lambda2*( Z3g*( 1 + 2*rho3^2) + 3*rho3*u3/sigma3) /
                      sqrt1r32^5
                      )
        hess <- array(NA, c( Nparam, Nparam))
        hess[igamma,igamma] <- l.gg
        hess[igamma,ibeta1] <- l.gb1; hess[ibeta1,igamma] <- t( l.gb1)
        hess[igamma,isigma1] <- l.gs2; hess[isigma1,igamma] <- t( l.gs2)
        hess[igamma,irho1] <- l.gr2; hess[irho1,igamma] <- t( l.gr2)
        hess[igamma,ibeta2] <- l.gb2; hess[ibeta2,igamma] <- t( l.gb2)
        hess[igamma,isigma2] <- l.gs3; hess[isigma2,igamma] <- t( l.gs3)
        hess[igamma,irho2] <- l.gr3; hess[irho2,igamma] <- t( l.gr3)
        hess[ibeta1,ibeta1] <- l.b1b1
        hess[ibeta1,isigma1] <- l.b1s2; hess[isigma1,ibeta1] <- t( l.b1s2)
        hess[ibeta1,irho1] <- l.b1r2; hess[irho1,ibeta1] <- t( l.b1r2)
        hess[ibeta1,ibeta2] <- 0; hess[ibeta2,ibeta1] <- 0
        hess[ibeta1,isigma2] <- 0; hess[isigma2,ibeta1] <- 0
        hess[ibeta1,irho2] <- 0; hess[irho2,ibeta1] <- 0
        hess[isigma1,isigma1] <- l.s2s2
        hess[isigma1,irho1] <- l.s2r2; hess[irho1,isigma1] <- l.s2r2
        hess[isigma1,ibeta2] <- 0; hess[ibeta2,isigma1] <- 0
        hess[isigma1,isigma2] <- 0; hess[isigma2,isigma1] <- 0
        hess[isigma1,irho2] <- 0; hess[irho2,isigma1] <- 0
        hess[irho1,irho1] <- l.r2r2
        hess[irho1,ibeta2] <- 0; hess[ibeta2,irho1] <- 0
        hess[irho1,isigma2] <- 0; hess[isigma2,irho1] <- 0
        hess[irho1,irho2] <- 0; hess[irho2,irho1] <- 0
        hess[ibeta2,ibeta2] <- l.b2b2
        hess[ibeta2,isigma2] <- l.b2s3; hess[isigma2,ibeta2] <- t( l.b2s3)
        hess[ibeta2,irho2] <- l.b2r3; hess[irho2,ibeta2] <- t( l.b2r3)
        hess[isigma2,isigma2] <- l.s3s3
        hess[isigma2,irho2] <- l.s3r3; hess[irho2,isigma2] <- t( l.s3r3)
        hess[irho2,irho2] <- l.r3r3
        hess
    }
    ## ---------------
    NXS <- ncol( XS)
    if(is.null(colnames(XS)))
        colnames(XS) <- rep("XS", NXS)
    NXO1 <- ncol( XO1)
    if(is.null(colnames(XO1)))
        colnames(XO1) <- rep("XO1", NXO1)
    NXO2 <- ncol( XO2)
    if(is.null(colnames(XO2)))
        colnames(XO2) <- rep("XO2", NXO2)
    Nparam <- NXS + NXO1 + NXO2 + 4
                                        # Total # of parameters
    Nobs <- length( YS)
    NO1 <- length( YS[YS==0])
    NO2 <- length( YS[YS==1])
    ## indices in for the parameter vector
    igamma <- 1:NXS
    ibeta1 <- seq(tail(igamma, 1)+1, length=NXO1)
    isigma1 <- tail(ibeta1, 1) + 1
    irho1 <- tail(isigma1, 1) + 1
    ibeta2 <- seq(tail(irho1, 1) + 1, length=NXO2)
    isigma2 <- tail(ibeta2, 1) + 1
    irho2 <- tail(isigma2, 1) + 1
    ## split data by selection
    YO1 <- YO1[YS==0]
    YO2 <- YO2[YS==1]
    XO1 <- XO1[YS==0,,drop=FALSE]
    XO2 <- XO2[YS==1,,drop=FALSE]
    XS1 <- XS[YS==0,,drop=FALSE]
    XS2 <- XS[YS==1,,drop=FALSE]
    if(print.level > 0) {
        cat( "Choice 2:", NO1, "times; choice 3:", NO2, "times\n")
    }
    if( print.level > 0) {
        cat( "Initial values:\n")
        cat("beta1\n")
        print(start[ibeta1])
        cat("sigma1\n")
        print(start[isigma1])
        cat("rho1\n")
        print(start[irho1])
        cat("beta2\n")
        print(start[ibeta2])
        cat("sigma2\n")
        print(start[isigma2])
        cat("rho2\n")
        print(start[irho2])
    }
    result <- maxLik(loglik, grad=gradlik, hess=hesslik, start=start,
                     print.level=print.level, ...)
#    compare.derivatives(loglik, gradlik, t0=start)
    class(result) <- c("tobit5", class(result))
    result
}
