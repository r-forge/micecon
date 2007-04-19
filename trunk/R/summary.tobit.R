summary.tobit <- function(object, ...) {
   ## object      object of class "selection" (only tobit2 and tobit5)
   ## ...         additional arguments for "summary.maxLik"
   ## 
   ## RESULTS:
   ## list of class "summary.selection" with following components:
   ## maximum    : function value at optimum
   ## estimateProbit : probit part of the estimate
   ## estimateEquation : equation -"-
   ## estimateEps      : sigma & rho part
   ## + additional components of "summary.maxLik"
   ##

   sl <- summary.maxLik(object, ...)
   sl$tobitType <- object$tobitType
   sl$estimateS <- sl$estimate[object$param$index$betaS,]

   if( object$tobitType == 2 ){
      sl$estimateO <- sl$estimate[object$param$index$betaO,]
      sl$estimateErr <- sl$estimate[c(object$param$index$sigma,
                                         object$param$index$rho),]
   } else if( object$tobitType == 5 ){
      sl$estimateO1 <- sl$estimate[object$param$index$betaO1,]
      sl$estimateO2 <- sl$estimate[object$param$index$betaO2,]
      sl$estimateErr <- sl$estimate[c(object$param$index$sigma1,
                                         object$param$index$sigma2,
                                         object$param$index$rho1,
                                         object$param$index$rho2),]
   }
   sl$param <- object$param
   return( sl )
}

print.summary.tobit <- function(x, ...) {

   cat("Tobit", x$tobitType, "model" )
   if( x$tobitType == 2 ) {
      cat( " (sample selection model)\n" )
   } else {
      cat( " (switching regression model)\n" )
   }
   cat( "Maximum Likelihood estimation\n" )
   cat(x$type, ", ", x$iterations, " iterations\n", sep="")
   cat("Return code ", x$code, ": ", x$message, "\n", sep="")
   if(!is.null(x$estimate)) {
      cat("Log-Likelihood:", loglikValue(x), "\n")
      cat( x$param$NObs, "observations" )
      if( x$tobitType == 2 ) {
         cat( " (", x$param$N0, " censored and ", x$param$N1, " observed)",
            sep = "" )
      } else {
         cat( " (", x$param$N1, " selection 1 and ",
            x$param$N2, " selection 2)", sep = "" )
      }
      cat( " and", x$param$NParam, "free parameters" )
      cat( " (df = ", x$param$df, ")\n", sep="")
      cat("Probit selection equation:\n")
      printCoefmat(x$estimateS, signif.legend=FALSE)
      if( x$tobitType == 2 ) {
         cat("Outcome equation:\n")
         printCoefmat(x$estimateO, signif.legend=FALSE)
      } else if( x$tobitType == 5 ) {
         cat("Outcome equation 1:\n")
         printCoefmat(x$estimateO1, signif.legend=FALSE)
         cat("Outcome equation 2:\n")
         printCoefmat(x$estimateO2, signif.legend=FALSE)
      }
      cat("Error terms:\n")
      printCoefmat(x$estimateErr)
   }
}
