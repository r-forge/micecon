cesEst <- function( yName, xNames, data, vrs = FALSE,
      method = "Nelder-Mead", ... ) {

   # y = gamma * ( alpha * x1^(-rho) + ( 1 - alpha ) * x2^(-rho) )^(-phi/rho)
   # s = 1 / ( 1 + rho )

   checkNames( c( yName, xNames ), names( data ) )

   if( length( xNames ) != 2 ) {
      stop( "currently, argument 'xNames' must contain exactly",
         " two variable names" )
   }

   # store the (matched) call
   matchedCall <- match.call()

   # prepare data for estimation
   estData <- data.frame( y = data[[ yName ]],
      x1 = data[[ xNames[ 1 ] ]], x2 = data[[ xNames[ 2 ] ]] )

   # start values
   startVal <- c( gamma = sqrt( mean( estData$y ) ),
      alpha = 0.5, rho = -0.5 )
   if( vrs ) {
      startVal <- c( startVal, phi = 1 )
   }

   # Estimation by the Kmenta approximation
   if( method == "Kmenta" ) {
      if( !vrs ) {
         warning( "allowing for variable returns to scale",
            " in the Kmanta approximation",
            " although argument 'vrs' is 'FALSE'." )
         matchedCall$vrs <- TRUE
      }
      result <- cesEstKmenta( yName = yName, xNames = xNames, data = data )
   } else if( method %in% c( "Nelder-Mead", "SANN" ) ) {
      result <- list()
      result$optim <- optim( par = startVal, fn = cesRss, data = estData,
         hessian = TRUE, method = method, ... )
      result$coefficients <- result$optim$par
   } else if( method %in% c( "BFGS", "CG", "L-BFGS-B" ) ) {
      result <- list()
      result$optim <- optim( par = startVal, fn = cesRss, gr = cesRssDeriv,
         data = estData, hessian = TRUE, method = method, ... )
      result$coefficients <- result$optim$par
   } else {
      stop( "argument 'method' must be either 'Nelder-Mead', 'BFGS',",
         " 'CG', 'L-BFGS-B', 'SANN', or 'Kmenta'" )
   }

   # covariance matrix of the estimated parameters
   if( is.null( result$vcov ) && !is.null( result$optim$hessian ) ) {
      if( det( result$optim$hessian ) >= .Machine$double.eps ) {
         result$vcov <- solve( result$optim$hessian )
      }
   }
   if( is.null( result$vcov ) ) {
      result$vcov <- matrix( NA, nrow = length( result$coefficients ),
         ncol = length( result$coefficients ) )
      rownames( result$vcov ) <- names( result$coefficients )
      colnames( result$vcov ) <- names( result$coefficients )
   }

   # return also the call
   result$call <- matchedCall

   # returned the method used for the estimation
   result$method <- method

   # fitted values
   result$fitted.values <- cesCalc( xNames = xNames, data = data,
      coef = result$coefficients )

   # nonlinear least squares
#    result$nls <- nls(
#       y ~ gamma * ( alpha * x1^rho + ( 1 - alpha ) * x2^rho )^(1/rho),
#       data = estData, start = result$startVal, trace = TRUE,
#       algorithm = "port", lower = c( -Inf, 0.01 , -Inf ),
#       upper = c( Inf, 0.99, Inf ) )

   class( result ) <- c( "cesEst", class( result ) )
   return( result )
}

