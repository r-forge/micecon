summary.tobit5 <- function(object, ...) {
   ## object      object of class "tobit5"
   ## ...         additional arguments for "summary.maxLik"
   ## 
   ## RESULTS:
   ## list of class "summary.tobit5" with following components:
   ## maximum    : function value at optimum
   ## estimateProbit : probit part of the estimate
   ## estimateEquation : equation -"-
   ## estimateEps      : sigma & rho part
   ## + additional components of "summary.maxLik"
   ## 
   sl <- summary.maxLik(object, ...)
   sl$estimateS <- sl$estimate[object$param$index$betaS,]
   sl$estimateO1 <- sl$estimate[object$param$index$betaO1,]
   sl$estimateO2 <- sl$estimate[object$param$index$betaO2,]
   sl$estimateErr <- sl$estimate[c(object$param$index$sigma1,
                                         object$param$index$sigma2,
                                         object$param$index$rho1,
                                         object$param$index$rho2),]
   sl$NObs <- object$param$NObs
   sl$N1   <- object$param$N1
   sl$N2   <- object$param$N2
   sl$NXS  <- object$param$NXS
   sl$NXO1 <- object$param$NXO1
   sl$NXO2 <- object$param$NXO2
   sl$df   <- object$param$df

   class(sl) <- c("summary.tobit5", class(sl))
   return( sl )
}

print.summary.tobit5 <- function(x, ...) {
   cat("--------------------------------------------\n")
   cat("Tobit 5 selection model/Maximum Likelihood estimation\n")
   cat(x$type, ", ", x$iterations, " iterations\n", sep="")
   cat("Return code ", x$code, ": ", x$message, "\n", sep="")
   if(!is.null(x$estimate)) {
      cat("Log-Likelihood:", loglikValue(x), "\n")
      cat(x$NObs, " observations (", x$N1, " selection 1 and ", x$N2, " selection 2) and ",
          x$NActivePar, " free parameters (df = ",
          x$NObs - x$NActivePar, ")\n", sep="")
      cat("Probit selection equation:\n")
      printCoefmat(x$estimateS, signif.legend=FALSE)
      cat("Outcome equation 1:\n")
      printCoefmat(x$estimateO1, signif.legend=FALSE)
      cat("Outcome equation 2:\n")
      printCoefmat(x$estimateO2, signif.legend=FALSE)
      cat("Error terms data:\n")
      printCoefmat(x$estimateErr)
   }
   cat("--------------------------------------------\n")
}
