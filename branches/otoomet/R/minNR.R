minNR <- function( f, theta0,
           names=NULL,
           print.level=0,
           tol=1e-6, gradtol=1e-6, steptol=1e-6,
           lambdatol=1e-6, iterlim=15, ...) {
   ### Newton-Raphson minimization. Note that Gauss-Newton is exactly as
   ### Newton-Raphson, except the hessian is approximated.
   ### Parameters:
   ### f - the function to be minimized
   ### theta0 - initial parameter vector
   ### names       - names of the parameters
   ### steptol     - minimum step size
   ### lambdatol   - max lowest eigenvalue when forcing pos. definite H
   ### ...         - extra arguments for f()
   ### The stopping criteria
   ### tol         - maximum allowed difference between sequential values
   ### gradtol     - maximum allowed norm of gradient vector
   ### iterlim     - maximum # of iterations
   ###
   ### results:
   ### samm        - antakse siis, kui tuleb sammu viga (tulemus
   ###             - 3). list, mis sisaldab:
   ###  teeta0     - parameetrid, millel viga tuli
   ###  f0         - funktsiooni väärtus nende parameetritega (koos
   ###               gradiendi ja hessi maatriksiga)
   ###  teeta1     - ruutpolünoomi järgi õige uus parameetri väärtus
   "minimization.message" <- function( code) {
      message <- switch(code,
         "1" = "gradient close to zero. May be a solution",
         "2" = paste("successive function values within tolerance",
            "limit.\n May be a solution"),
         "3" = paste("Last step could not find a value below the",
            "current.\nMay be near a solution"),
         "4" = "Iteration limit exceeded.",
         "100" = "Initial value out of range.",
         paste("Code", code))
      return(message)
   }
   min.eigen <- function( M) {
      ## return minimal eigenvalue of (symmetric) matrix
      eigen( M, symmetric=TRUE, only.values=TRUE)$values[K]
      ## L - eigenvalues in decreasing order
   }

   maximisation.type <- "Newton-Raphson minimization"
   K <- length( theta0)
   I <- diag( seq( 1, 1, length=K))
                                       # I is unit matrix
   theta1 <- theta0
   iter <- 0
   f1 <- f( theta1, ...)
   if(is.na( f1)) {
      result <- list(
                     code=100, message=minimization.message(100),
                     iterations=iter,
                     type=maximisation.type)
      class(result) <- "maximisation"
      return(result)
   }
   G <- attr( f1, "gradient")
   if( is.null( G)) {
      stop( "Jama: puudub gradient\n")
   }
   if( is.null( attr( f1, "hessian"))) {
      stop( "Jama: puudub hessi maatriks\n")
   }
   repeat {
      iter <- iter + 1
      lambda <- NA
      theta <- theta1
      f0 <- f1
      H0 <- attr( f0, "hessian")
                                        # Now enforce H to be positive definite
      H <- H0
      step <- 1
      ## Note: if we are increasing lambda, we should increase max
      ## allowed step too
      if( min.eigen( H) <= lambdatol) {
         ## If smallest eigenvalue is very small, matrix is near singular
         H <- H0 + I
         ## try lambda=1 and upward first
         lambda <- 1
         if( min.eigen( H) <= lambdatol) {
            while( min.eigen( H) <= lambdatol) {
               lambda <- 2*lambda
               H <- H0 + lambda*I
            }
         } else {
            ## try to find a smaller lambda
            while( min.eigen( H) > lambdatol) {
               lambda <- lambda/2
               step <- step*2
               H <- H0 + lambda*I
            }
            lambda <- lambda*2
            H <- H0 + lambda*I
         }
      }
      ## Now H should be a relatively small positive definite matrix,
      ## somehow connected with hessian
      iH <- solve( H)
                                        # solve() inverts matrix
      amount <- iH %*% G
      theta1 <- theta - step*amount
                                        # Note: ,-' means we are
                                        # climbing downhill
      f1 <- f(theta1, ...)
      while( is.na( f1) || ( ( f1 >= f0) && ( step > steptol))) {
         # We end up in a NA or a higher value. try smaller step
         step <- step/2
         theta1 <- theta - step*amount
         f1 <- f(theta1, ...)
      }
      G <- attr( f1, "gradient")
      if( print.level > 1) {
         cat( "-----Iteration", iter, "-----\n")
         cat( "lambda ", lambda, " step", step, " fcn value:",
            as.vector(f1), "\n")
         a <- cbind(G, amount, theta1)
         dimnames(a) <- list(names, c("gradient", "amount", "new param"))
         print(a)
         if( print.level > 2) {
            H <- attr( f1, "hessian")
            cat( "Hessian: det=", det( H), "\n")
            print( H)
         }
      }
      if( step < steptol) {
         code <- 3; break
      }
      if( iter > iterlim) {
         code <- 4; break
      }
      if( sqrt( t( G)%*%G) < gradtol) {
         code <-1; break
      }
      if( f0 - f1 < tol) {
         code <- 2; break
      }
   }
   if( print.level > 0) {
      cat( "--------------\n")
      cat( minimization.message( code))
      cat( iter, " iterations\n")
      cat( "estimate:", theta1, "\n")
      cat( "Function value:", f1, "\n")
   }
   if( code == 3) {
      samm <- list(theta0=theta, f0=f0, theta1=theta1)
   } else {
      samm <- NULL
   }
   H <- attr( f1, "hessian")

   result <-list(
                  maximum=as.vector( f1),
                  estimate=theta1,
                  gradient=-G,
                  hessian=-H,
                  code=code,
                  message=minimization.message( code),
                  last.step=samm,
                                        # only when could not find a
                                        # lower point
                  iterations=iter,
                  type=maximisation.type)
   class(result) <- "maximisation"
   invisible(result)
}

