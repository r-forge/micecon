## Return Hessian of an object

Hessian <- function(x, ...)
    UseMethod("Hessian")

Hessian.default <- function(x, ...)
    x$Hessian
