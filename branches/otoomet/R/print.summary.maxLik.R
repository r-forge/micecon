print.summary.maxLik <- function( x, ... ) {
   cat("--------------------------------------------\n")
   cat("Maximum Likelihood estimation\n")
   cat(x$type, ", ", x$iterations, " iterations\n", sep="")
   cat("Return code ", x$code, ": ", x$message, "\n", sep="")
   if(!is.null(x$estimate)) {
      cat("Log-Likelihood:", x$loglik, "\n")
      cat(x$NActivePar, " free parameters\n")
      cat("Estimates:")
      print(x$estimate)
      if(!is.null(x$Hessian)) {
         cat("Hessian:\n")
         print(x$Hessian)
      }
   }
   cat("--------------------------------------------\n")
}

