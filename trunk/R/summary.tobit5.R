
summary.tobit5 <- function(object, ...) {
   sl <- summary.maxLik(object, ...)
   s <- c(sl,
          estimateS=list(sl$estimate[object$param$index$betaS,]),
          estimateO1=list(sl$estimate[object$param$index$betaO1,]),
          estimateO2=list(sl$estimate[object$param$index$betaO2,]),
          estimateErr=list(sl$estimate[c(object$param$index$sigma1,
                                         object$param$index$rho1,
                                         object$param$index$sigma2,
                                         object$param$index$rho2),]),
          NObs=object$param$NObs, NActivePar=object$param$NActivePar,
          N1=object$param$N1, N2=object$param$N2,
          NXS=object$param$NXS, NXO1=object$param$NXO, NXO2=object$param$NXO2,
          df=object$df
          )
   class(s) <- c("summary.tobit5", class(sl))
   s
}

print.summary.tobit5 <- function(x, ...) {
   cat("--------------------------------------------\n")
   cat("Tobit 5 selection model/Maximum Likelihood estimation\n")
   cat(x$type, ", ", x$iterations, " iterations\n", sep="")
   cat("Return code ", x$code, ": ", x$message, "\n", sep="")
   if(!is.null(x$estimate)) {
      cat("Log-Likelihood:", loglikValue(x), "\n")
      cat(x$NObs, " observations (", x$N1, " selection 1 ", x$N2, " selection 2) and ",
          x$NActivePar, " free parameters (df =",
          x$NObs - x$NActivePar, ")\n", sep="")
      cat("\nProbit selection equation:\n")
      print(x$estimateS)
      cat("\nOutcome 1:\n")
      print(x$estimateO1)
      cat("\nOutcome 2:\n")
      print(x$estimateO2)
      cat("\nError terms data:\n")
      print(x$estimateErr)
   }
   cat("--------------------------------------------\n")
}
