## Return the #of parameters of model
NParam <- function(x, ...)
    UseMethod("NParam")

NParam.default <- function(x, ...)
    x$param$NParam

NParam.lm <- function(x, ...)
    length(coefficients(x))
