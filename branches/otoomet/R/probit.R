### -*- ess-indent-level: 3 -*-
probit <- function( formula, b0=NULL, data=sys.frame(sys.parent()),
                   ...) {
   ## Probit binary choice model
   ## formula: model formula, response must be either a logical or numeric vector containing only 0-s and
   ##          1-s 
   ## b0:      initial value of the parameters
   ## ...      further arguments for the maxLik algorithm
   ##
   ## return: a list with following components.
   ##  $results: maximisation results
   ##  $LRT:     list with two components:
   ##            LRT: likelihood ration test H0 - none of the variables significant
   ##            df:  corresponding degrees of freedom
   ##
   loglik <- function( beta) {
      xb0 <- x0 %*% beta
      xb1 <- x1 %*% beta
      loglik <- sum(pnorm( xb0, log=TRUE, lower.tail=FALSE)) + sum(pnorm( xb1, log.p=TRUE))
   }
   gradlik <- function(beta) {
      xb0 <- x0 %*% beta
      xb1 <- x1 %*% beta
      gradlik <- - t(x0) %*% (dnorm(xb0)/pnorm(xb0, lower.tail=FALSE)) +
          t(x1) %*% (dnorm(xb1)/pnorm(xb1))
   }
   hesslik <- function(beta) {
      xb0 <- as.vector( x0 %*% beta)
      xb1 <- as.vector( x1 %*% beta)
      F0 <- pnorm( xb0)
      F1 <- pnorm( xb1)
      f0 <- dnorm( xb0)
      f1 <- dnorm( xb1)
      yF0 <- pnorm( xb0, lower.tail=FALSE)
      loglik2 <- t( x1) %*% ( x1 * ( -f1*xb1*F1 - f1*f1)/F1/F1) -
         t( x0) %*% ( x0 * ( -f0*xb0*yF0 + f0*f0)/yF0/yF0)
                    # note that df/db' = -f (x'b) x'
   }
   x <- model.matrix(formula, data=data)
   y <- as.numeric(model.response( model.frame( formula, data=data)))
   NParam <- ncol( x)
   NObs <- length( y)
   N1 <- length( y[y != 0])
   N0 <- NObs - N1
   if(N0 == 0 | N1 == 0) {
      stop("only zero or one observations")
   }
   x0 <- x[y==0,]
   x1 <- x[y==1,]
   if(is.null(b0)) {
      b0 <- rep( 0, NParam)
   }
   if(is.null(names(b0))) {
      names(b0) <- dimnames(x)[[2]]
   }
   ## Main estimation
   estimation <- maxLik(loglik, gradlik, hesslik, b0,
                        method="Newton-Raphson", ...)
   ## compare.derivatives(gradlik, hesslik, t0=b0)
                                        #
   ## Likelihood ratio test: H0 -- all the coefficients, except intercept
   ## are zeros.  ML estimate for this model is qnorm(N1/NObs)
   ll.bar <- loglik(c(qnorm(N1/NObs), rep(0, NParam-1)))
   LRT <- 2*(estimation$maximum - ll.bar)
                                        #
   result <- c(estimation,
               LRT=list(list(LRT=LRT, df=NParam-1)),
               NParam=NParam,
               NObs=NObs,
               N1=N1,
               N0=N0,
               df=NObs - NParam)
                                        # there are df-1 constraints
   class(result) <- c("probit", class(estimation))
   result
}
