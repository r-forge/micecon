tobit2 <- function(selection, reg, data=sys.frame(sys.parent()),
                   init=NULL, print.level=0, ...) {
   ### The model (Amemiya 1985):
   ### The latent variables are:
   ### y1* = z'gamma + u1
   ### y2* = x'beta + u2
   ### The observables are:
   ###      / 1  if  y1* > 0
   ### y1 = \ 0  if  y1* <= 0
   ###      / 0    if  y1 = 0
   ### y2 = \ y2*  if  y1 = 1
   ###
   ### PARAMETERS:
   ### selection: formula for the selection equation (y1)
   ### reg        main model
   ### data       dataset
   ### init       initial value of coefficients.  The order is as follows:
   ###            init = (gamma', beta', sigma, rho)'
   ### print.level: 0 - nothing printed, as larger, as more information printed
   ###  ...         additional parameters form maximisation
   ###
   ### RESULTS:
   ### a list with the following components:
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
      u2 <- Y2 - X2.b
      r <- sqrt( 1 - rho^2)
      B <- (Z2.g + rho/sigma*u2)/r
      l1 <- sum(pnorm(-Z1.g, log.p=TRUE))
      # loglik related with unobserved cases
      l2 <- -N2/2*log(2*pi) - N2*log(sigma) +
         sum(-0.5*(u2/sigma)^2 + pnorm(B, log.p=TRUE))
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
      u2 <- Y2 - X2.b
      r <- sqrt( 1 - rho^2)
      B <- (Z2.g + rho/sigma*u2)/r
      lambdaB <- dnorm(B)/pnorm(B)
      gradient <- numeric(NParam)
      gradient[igamma] <- t(Z1) %*% (-dnorm(-Z1.g)/pnorm(-Z1.g)) +
         (t(Z2) %*% lambdaB)/r
      gradient[ibeta] <- t(X2) %*% (u2/sigma^2 - lambdaB*rho/sigma/r)
      gradient[isigma] <- sum(u2^2/sigma^3 - 1/sigma - lambdaB*rho*u2/sigma^2/r)
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
      u2 <- Y2 - X2.b
      r <- sqrt( 1 - rho^2)
      B <- (Z2.g + rho/sigma*u2)/r
      lambdaB <- dnorm(B)/pnorm(B)
      fZ1.g <- dnorm(-Z1.g)
      FZ1.g <- pnorm(-Z1.g)
      C <- -(pnorm(B)*dnorm(B)*B + dnorm(B)^2)/pnorm(B)^2
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
                                   1/sigma^2
                                   -3*u2*u2/sigma^4
                                   +2*lambdaB* u2/r *rho/sigma^3
                                   +rho^2/sigma^4 *u2*u2/r^2 *C)
      hess[isigma,irho] <- hess[irho,isigma] <-
         -sum((C*rho*(u2/sigma + rho*Z2.g)/r + lambdaB)*
         u2/sigma^2)/r^3
      hess[irho,irho] <-
         sum(C*((u2/sigma + rho*Z2.g)/r^3)^2 +
         lambdaB*(Z2.g*(1 + 2*rho^2) + 3*rho*u2/sigma) / r^5 )
      return( hess )
   }
   ## --- the main program ---
   Z <- model.matrix(selection, data=data)
   y1 <- model.response(model.frame(selection, data=data))
   NZ <- ncol( Z)
   X <- model.matrix(reg, data=data)
   y2 <- model.response(model.frame(reg, data=data))
   NX <- ncol( X)
   NParam <- NZ + NX + 2
                                        # Total # of parameters
   NObs <- length( y1)
   N1 <- length( y1[y1==0])
   N2 <- length( y1[y1==1])
   ## indices in for the parameter vector
   igamma <- 1:NZ
   ibeta <- max(igamma) + seq(length=NX)
   isigma <- max(ibeta) + 1
   irho <- max(isigma) + 1
   ## divide data by choices
   Y2 <- y2[y1==1]
   X2 <- X[y1==1,,drop=FALSE]
   Z1 <- Z[y1==0,,drop=FALSE]
   Z2 <- Z[y1==1,,drop=FALSE]
   if(print.level > 0) {
      cat( "Not observed:", N1, "; observed:", N2,"\n", sep="")
   }
   if(is.null(init)) {
      init <- numeric(NParam)
      ## initial values.  First probit part w/probit model
      probit <- probit(selection, data=data, print.level=print.level-2)
      if( print.level > 0) {
         cat("The probit part of the model:\n")
         print(summary(probit))
      }
      gamma <- probit$results$estimate
      init[igamma] <- gamma
      ## regression parts by Heckman's two-step method
      lambda <- dnorm( Z2%*%gamma)/pnorm( Z2%*%gamma)
      delta <- mean( lambda^2 + Z2%*%gamma *lambda)
      twoStep <- summary( lm( Y2 ~ X2 + lambda))
      init[ibeta] <- twoStep$coefficients[1:NX,1]
      se <- twoStep$sigma
      beta.l <- twoStep$coefficients["lambda",1]
      sigma <- sqrt(se^2 + beta.l^2*delta)
      init[isigma] <- sigma
      rho <- beta.l/sigma
      if( rho <= -0.99) rho <- -0.99
      if( rho >= 0.99) rho <- 0.99
      init[irho] <- rho
   }
   if(is.null(names(init))) {
      names(init) <- c(colnames(Z), colnames(X), "sigma", "rho")
   }
   if( print.level > 0) {
      cat( "Initial values (Heckman two-step estimate):\n")
      print(init)
   }
   results <- maxNR(loglik, gradlik, theta=init,
                     print.level=print.level - 1, ...)
   # compare.derivatives(gradlik, hesslik, t0=init, ...)
   tobit2 <- list(probit=probit,
                   results=results)
   class(tobit2) <- "tobit2"
   return( tobit2 )
}
