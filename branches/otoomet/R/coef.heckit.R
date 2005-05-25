coef.heckit <- function(x) {
   b <- c(coef(x$probit), coef(x$lm), x$sigma, x$rho)
   N <- length(b)
   names(b)[(N-1):N] <- c("sigma", "rho")
   b
}
