print.summary.maximisation <- function( x, ... ) {
   summary <- x
   cat("--------------------------------------------\n")
   cat(summary$type, "\n")
   cat("Number of iterations:", summary$iterations, "\n")
   cat("Return code:", summary$code, "\n")
   cat(summary$message, "\n")
   if(!is.null(summary$unsucc.step)) {
      cat("Last (unsuccessful) step: function value", summary$unsucc.step$value,
         "\n")
      print(summary$unsucc.step$parameters)
   }
   if(!is.null(summary$estimate)) {
      cat("Function value:", summary$estimate$value, "\n")
      cat("Estimates:")
      print(summary$estimate)
      if(!is.null(summary$Hessian)) {
         cat("Hessian:\n")
         print(summary$Hessian)
      }
   }
   cat("--------------------------------------------\n")
}

