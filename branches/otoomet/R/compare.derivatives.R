compare.derivatives <- function(f, grad, hess=NULL, t0, eps=1e-6, ...) {
   ### t0 - initial parameter vector
   cat("-------- compare derivatives -------- \n")
   a <- f(t0, ...)
   analytic <- grad(t0, ...)
   if(is.null(dim(analytic))) {
      cat("Note: analytic gradient is vector.  Transforming into a matrix form\n")
      analytic <- matrix(analytic, 1, length(analytic))
   }
   cat("Dim of analytic gradient:", dim(analytic), "\n")
   numeric <- numeric.gradient(f, t0, eps, ...)
   cat("       numeric          :", dim(numeric), "\n")
   a <- rbind(t0, analytic, numeric, (analytic - numeric)/analytic)
   dimnames(a) <- list(c("theta 0", "analytic", "numeric", "rel.diff"),
                   param=NULL)
   print(a)
   if(!is.null(hess)) {
      analytic <- hess(t0, ...)
      numeric <- numeric.gradient(grad, t0, eps, ...)
      print((analytic - numeric)/analytic)
   }
   cat("-------- END of compare derivatives -------- \n")
}
