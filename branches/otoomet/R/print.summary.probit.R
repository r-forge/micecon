print.summary.probit <- function(x, ...) {
   cat("--------------------------------------------\n")
   cat("Probit model\n")
   cat("Number of iterations:", x$opt$iterations, "\n")
   cat("Return code:", x$opt$code, "\n")
   cat(x$opt$message, "\n")
   if(!is.null(x$opt$estimate)) {
      cat("Log-likelihood:", x$opt$estimate$value, "\n")
      cat("Estimates:")
      print(x$opt$estimate$results)
   }
   cat("LRT H0: no dependence on covariates\n")
   cat("chi2(", x$LRT$df, ") =", x$LRT$LRT, " (p=", x$LRT$pchi2, ")\n")
   cat("--------------------------------------------\n")
}
