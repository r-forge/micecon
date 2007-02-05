coef.heckit <- function(object, ...) {
   b <- c(coef(object$probit), object$step2coef, object$sigma, object$rho)
   N <- length(b)
   names(b)[(N-1):N] <- c("sigma", "rho")
   b
}
