tobit2 <- function(selection, formula, method="ml",
                   data=sys.frame(sys.parent()),
                   b0=NULL, print.level=0, ...) {
   ## The model (Amemiya 1985):
   ## The latent variables are:
   ## y1* = z'gamma + u1
   ## y2* = x'beta + u2
   ## The observables are:
   ##      / 1  if  y1* > 0
   ## y1 = \ 0  if  y1* <= 0
   ##      / 0    if  y1 = 0
   ## y2 = \ y2*  if  y1 = 1
   ##
   ## PARAMETERS:
   ## selection: formula for the selection equation (y1)
   ## formula    main model
   ## method     maximum likelihood (ml) or two-step (2step)
   ## data       dataset
   ## b0       initial value of coefficients.  The order is as follows:
   ##            b0 = (gamma', beta', sigma, rho)'
   ## print.level: 0 - nothing printed, as larger, as more information printed
   ##  ...         additional parameters form maximisation
   ##
   ## RESULTS:
   ## a list with the following components:
   ##
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
   ## heckit -- two-step method
   heckit <- function( selection, formula, data ) {
      result <- list()
      result$probit <- probit( selection, data=data, x=TRUE)
      
      data$probitLambda <- dnorm( result$probit$linear.predictors ) /
          pnorm( result$probit$linear.predictors )

      data$probitDelta <- data$probitLambda * ( data$probitLambda +
                                               result$probit$linear.predictors )

      step2formula <- as.formula( paste( formula[ 2 ], "~", formula[ 3 ],
                                        "+ probitLambda" ) )

      result$lm <- lm( step2formula, data, data$probitdummy == 1 )

      result$sigma <- as.numeric( sqrt( crossprod( residuals( result$lm ) ) /
                                       sum( data$probitdummy == 1 ) +
                                       mean( data$probitDelta[ data$probitdummy == 1 ] ) *
                                       coefficients( result$lm )[ "probitLambda" ]^2 ) )

      result$rho <- coefficients( result$lm )[ "probitLambda" ] / result$sigma
      result$probitLambda <- data$probitLambda
      result$probitDelta  <- data$probitDelta

                                        # the foolowing variables are named according to Greene (2003), p. 785
      xMat <- model.matrix( result$lm )
      wMat <- model.matrix( result$probit )[ data$probitdummy == 1, ]
      fMat <- t( xMat ) %*% diag( result$probitDelta[
                                                     data$probitdummy == 1 ] ) %*% wMat
      qMat <- result$rho^2 * ( fMat %*% vcov( result$probit )%*% t( fMat ) )
      result$vcov <- result$sigma^2 * solve( crossprod( xMat ) ) %*%
          ( t( xMat ) %*% diag( 1 - result$rho^2 *
                               result$probitDelta[ data$probitdummy == 1 ] ) %*%
           xMat + qMat ) %*% solve( crossprod( xMat ) )
      result$coef <- matrix( NA, nrow = length( coef( result$lm ) ), ncol = 4 )
      rownames( result$coef ) <- names( coef( result$lm ) )
      colnames( result$coef ) <- c( "Estimate", "Std. Error", "t value",
                                   "Pr(>|t|)" )
      result$coef[ , 1 ] <- coef( result$lm )
      result$coef[ , 2 ] <- sqrt( diag( result$vcov ) )
      result$coef[ , 3 ] <- result$coef[ , 1 ] / result$coef[ , 2 ]
      result$coef[ , 4 ] <- 2 * ( 1 - pt( abs( result$coef[ , 3 ] ),
                                         result$lm$df ) )

      class( result ) <- "heckit"
      return( result )
   }
   ## --- the main program ---
   ## First the consistency checks
   if( class( formula ) != "formula" ) {
      stop( "argument 'formula' must be a formula" )
   }
   if( length( formula ) != 3 ) {
      stop( "argument 'formula' must be a 2-sided formula" )
   }
   if( "probit" %in% substr( all.vars( formula ), 1, 6 ) ) {
      stop( paste( "argument 'formula' may not include variable names",
                  "starting with 'probit'" ) )
   }
   if( class( selection ) != "formula" ) {
      stop( "argument 'selection' must be a formula" )
   }
   if( length( selection ) != 3 ) {
      stop( "argument 'selection' must be a 2-sided formula" )
   }
   if( "probit" %in% substr( all.vars( selection ), 1, 6 ) ) {
      stop( paste( "argument 'selection' may not include a variable",
                  "names starting with 'probit'" ) )
   }
   data$probitdummy <- model.frame( selection, data = data )[ , 1 ]
   test <- levels( as.factor( as.numeric( data$probitdummy ) ) )
   if( length( test ) != 2 ) {
      stop( paste( "The left hand side of 'selection' may only contain",
                  "1 and 0 or TRUE and FALSE" ) )
   }
   if( !all.equal( test, c( "0", "1" ) ) ) {
      stop( paste( "The left hand side of 'selection' may only contain",
                  "1 and 0 or TRUE and FALSE" ) )
   }
   ## now check whether two-step method is needed: either for final estimate or initial parameters
   if(method == "2step" | is.null(b0)) {
      twoStep <- heckit(selection, formula, data)
      if(method == "2step") {
         return(twoStep)
      }
      b0 <- coefficients(twoStep)
      if(print.level > 0) {
         cat("two-step estimate:\n")
         print(summary(twoStep))
      }
   }
   Z <- model.matrix(selection, data=data)
   y1 <- model.response(model.frame(selection, data=data))
   NZ <- ncol( Z)
   X <- model.matrix(formula, data=data)
   y2 <- model.response(model.frame(formula, data=data))
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
   results <- maxLik(loglik, gradlik, hesslik,
                     theta=b0,
                     print.level=print.level - 1, ...)
   # compare.derivatives(gradlik, hesslik, t0=b0, ...)
   results$twoStep <- twoStep
   results$NZ <- NZ
   results$NX <- NX
   results$N1 <- N1
   results$N2 <- N2
   class(results) <- c("tobit2", class(results))
   return(results)
}
