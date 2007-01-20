summary.tobit2 <- function(object, ...) {
   ## object      object of class "tobit2"
   ## ...         additional arguments for "summary.maxLik"
   ## 
   ## RESULTS:
   ## list of class "summary.tobit2" with following components:
   ## maximum    : function value at optimum
   ## estimateProbit : probit part of the estimate
   ## estimateEquation : equation -"-
   ## estimateEps      : sigma & rho part
   ## + additional components of "summary.maxLik"
   ## 
   sl <- summary.maxLik(object, ...)
   s <- c(sl,
          estimateS=list(sl$estimate[object$param$index$betaS,]),
          estimateO=list(sl$estimate[object$param$index$betaO,]),
          estimateErr=list(sl$estimate[c(object$param$index$sigma,
                                         object$param$index$rho),]),
          NObs=object$param$NObs, NActivePar=object$param$NActivePar,
          N0=object$param$N0, N1=object$param$N1,
          NXS=object$param$NXS, NXO=object$param$NXO, df=object$param$df
          )
   class(s) <- c("summary.tobit2", class(sl))
   s
}

print.summary.tobit2 <- function(x, ...) {
   cat("--------------------------------------------\n")
   cat("Tobit 2 selection model/Maximum Likelihood estimation\n")
   cat(x$type, ", ", x$iterations, " iterations\n", sep="")
   cat("Return code ", x$code, ": ", x$message, "\n", sep="")
   if(!is.null(x$estimate)) {
      cat("Log-Likelihood:", loglikValue(x), "\n")
      cat(x$NObs, " observations (", x$N0, " censored and ", x$N1, " observed) and ",
          x$NActivePar, " free parameters (df =",
          x$NObs - x$NActivePar, ")\n", sep="")
      cat("\nProbit selection equation:\n")
      print(x$estimateS)
      cat("\nOLS equation:\n")
      print(x$estimateO)
      cat("\nError terms data:\n")
      print(x$estimateErr)
   }
   cat("--------------------------------------------\n")
}
