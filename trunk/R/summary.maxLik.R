### summary methods for 'maxLik' class, i.e. 'raw' ML estimation results.  Class 'MLEstimate' is for
### complete models.

print.summary.maxLik <- function( x, ... ) {
   cat("--------------------------------------------\n")
   cat("Maximum Likelihood estimation\n")
   cat(x$type, ", ", x$iterations, " iterations\n", sep="")
   cat("Return code ", returnCode(x), ": ", x$message, "\n", sep="")
   if(!is.null(x$estimate)) {
      cat("Log-Likelihood:", x$loglik, "\n")
      cat(x$NActivePar, " free parameters\n")
      cat("Estimates:\n")
      print(x$estimate)
      if(!is.null(Hessian(x))) {
         cat("Hessian:\n")
         print(Hessian(x))
      }
   }
   cat("--------------------------------------------\n")
}

summary.maxLik <- function( object, hessian=FALSE, ... ) {
   ## object      object of class "maxLik"
   ## hessian     logical, whether to include hessian in summary
   ## 
   ## RESULTS:
   ## list of class "summary.maxLik" with following components:
   ## maximum    : function value at optimum
   ## estimate   : estimated parameter values at optimum
   ## gradient   :           gradient at optimum
   ## hessian    :           hessian
   ## code       : code of convergence
   ## message    : message, description of the code
   ## iterations : number of iterations
   ## type       : type of optimisation
   ##
   result <- object$maximisation
   NParam <- length(coef <- coefficients(object))
   if(!is.null(object$activePar)) {
      activePar <- object$activePar
   } else {
      activePar <- rep(TRUE, NParam)
   }
   if(returnCode(object) < 100) {
      if(min(abs(eigen(Hessian(object)[activePar,activePar],
                       symmetric=TRUE, only.values=TRUE)$values)) > 1e-6) {
         varcovar <- matrix(0, NParam, NParam)
         varcovar[activePar,activePar] <-
             solve(-Hessian(object)[activePar,activePar])
         hdiag <- diag(varcovar)
         if(any(hdiag < 0)) {
            warning("Diagonal of variance-covariance matrix not positive!\n")
         }
         stdd <- sqrt(hdiag)
         t <- coef/stdd
         p <- 2*pnorm( -abs( t))
      } else {
         stdd <- 0
         t <- 0
         p <- 0
      }
      results <- cbind(coef, stdd, t, "P(|b| > t)"=p)
      Hess <- NULL
      if(hessian) {
         Hess <- Hessian(object)
      }
   } else {
      results <- NULL
      Hess <- NULL
   }
   summary <- list(type=object$type,
                   iterations=object$iterations,
                   code=returnCode(object),
                   message=object$message,
                   loglik=object$maximum,
                   estimate=results,
                   Hessian=Hess,
                   activePar=object$activePar,
                   NActivePar=sum(object$activePar))
   class(summary) <- "summary.maxLik"
   summary
}
