numeric.nHessian <- function( f, t00, eps=1e-6, ...) {
### Numeric Hessian without gradient
    f00 <- f( t00, ...)
    eps2 <- eps*eps
    N <- length( t00)
    H <- matrix(NA, N, N)
    for( i in 1:N) {
      for( j in 1:N) {
        t01 <- t00
        t10 <- t00
        t11 <- t00
                                        # initial point
        t01[i] <- t01[i] + eps
        t10[j] <- t10[j] + eps
        t11[i] <- t11[i] + eps
        t11[j] <- t11[j] + eps
        f01 <- f( t01, ...)
        f10 <- f( t10, ...)
        f11 <- f( t11, ...)
        H[i,j] <- ( f11 - f01 - f10 + f00)/eps2
      }
    }
    H
  }
