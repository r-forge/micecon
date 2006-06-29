## Return #of observations for models

NObs <- function(x, ...)
    ## Number of observations for statistical models
    UseMethod("NObs")

NObs.lm <- function(x, ...)
    nrow(x$qr$qr)

NObs.probit <- function(x)
    x$param$NObs
