numericHessian <- function(f, grad=NULL, t0, eps=1e-6, ...) {
   a <- f(t0, ...)
   if(is.null(grad)) {
      numeric.nHessian(f, t0, eps, ...)
   } else {
      numericGradient(grad, t0, eps, ...)
   }
}
