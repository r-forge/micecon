print.summary.tobit5 <- function(x, ...) {
    cat("--------------------------------------------\n")
    cat("Tobit 5 model\n")
    cat("Number of iterations:", x$opt$iterations, "\n")
    cat("Return code:", x$opt$code, "(", x$opt$message, ")\n")
    if(!is.null(x$opt$estimate)) {
        cat("Log-likelihood:", x$opt$estimate$value, "\n")
        cat("Estimates:")
        print(x$opt$estimate$results)
    }
    cat("--------------------------------------------\n")
}
