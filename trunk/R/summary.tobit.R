summary.tobit <- function(object, ...) {
   ## object      object of class "selection" (only tobit2 and tobit5)
   ## ...         additional arguments for "summary.maxLik"

   s <- summary.maxLik(object, ...)
   s$tobitType <- object$tobitType
   s$param <- object$param
   return( s )
}

print.summary.tobit <- function(x,
                                 part="full",
                                 ... ) {

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
      if(part == "full") {
         cat("Probit selection equation:\n")
         printCoefmat( x$estimate[ x$param$index$betaS, ],
            signif.legend = FALSE )
      }
      if( x$tobitType == 2 ) {
         cat("Outcome equation:\n")
         printCoefmat( x$estimate[ x$param$index$betaO, ],
            signif.legend = ( part != "full" ) )
      } else if( x$tobitType == 5 ) {
         cat("Outcome equation 1:\n")
         printCoefmat( x$estimate[ x$param$index$betaO1, ],
            signif.legend = FALSE )
         cat("Outcome equation 2:\n")
         printCoefmat( x$estimate[ x$param$index$betaO2, ],
            signif.legend = ( part != "full" ) )
      }
      if(part=="full") {
         if( x$tobitType == 2 ) {
            i <- c(x$param$index$Mills, x$param$index$sigma, x$param$index$rho)
         } else if( x$tobitType == 5 ) {
            i <- c( x$param$index$Mills1, x$param$index$Mills2,
               x$param$index$sigma1, x$param$index$sigma2,
               x$param$index$rho1, x$param$index$rho2 )
         }
         cat("Error terms:\n")
         printCoefmat(x$estimate[i,])
      }
   }
}
