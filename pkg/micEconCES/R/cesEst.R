cesEst <- function( yName, xNames, data, vrs = FALSE,
      method = "Nelder-Mead",
      startVal = c( sqrt( mean( data[[ yName ]] ) ), 0.5, -0.5, 1 ),
      ... ) {

   # y = gamma * ( delta * x1^(-rho) + ( 1 - delta ) * x2^(-rho) )^(-phi/rho)
   # s = 1 / ( 1 + rho )

   checkNames( c( yName, xNames ), names( data ) )

   if( length( xNames ) != 2 ) {
      stop( "currently, argument 'xNames' must contain exactly",
         " two variable names" )
   }

   # start values
   if( length( startVal ) != 4 && vrs ) {
      stop( "a CES function with 2 explanatory variables and VRS has 4",
         " parameters but you provided ", length( startVal ),
         " starting values" )
   } else if( ! length( startVal ) %in% c( 3, 4 ) && !vrs ) {
      stop( "a CES function with 2 explanatory variables and CRS has 3",
         " parameters but you provided ", length( startVal ),
         " starting values" )
   }
   names( startVal ) <- c( "gamma", "delta", "rho", "phi" )[
      1:length( startVal ) ]
   startVal <- startVal[ 1:( 3 + vrs ) ]

   # store the (matched) call
   matchedCall <- match.call()

   # prepare data for estimation
   estData <- data.frame( y = data[[ yName ]],
      x1 = data[[ xNames[ 1 ] ]], x2 = data[[ xNames[ 2 ] ]] )

   # prepare list that will be returned
   result <- list()

   # Estimation by the Kmenta approximation
   if( method == "Kmenta" ) {
      result <- cesEstKmenta( yName = yName, xNames = xNames, data = data,
         vrs = vrs )
   } else if( method %in% c( "Nelder-Mead", "SANN", "BFGS", "CG", "L-BFGS-B" ) ) {
      if( method %in% c( "Nelder-Mead", "SANN" ) ) {
         result$optim <- optim( par = startVal, fn = cesRss, data = estData,
            method = method, ... )
      } else {
         result$optim <- optim( par = startVal, fn = cesRss, gr = cesRssDeriv,
            data = estData, method = method, ... )
      }
      result$coefficients <- result$optim$par
      result$iter <- min( result$optim$counts, na.rm = TRUE )
      result$convergence <- result$optim$convergence == 0
      result$message <- result$optim$message
   } else if( method == "LM" ) {
      # residual function
      residFun <- function( par, data2 ) {
         result <- data2$y - cesCalc( xNames = c( "x1", "x2" ),
            data = data2, coef = par )
         return( result )
      }

      # jacobian function
      jac <- function( par, data3 ) {
         return( -c( cesDerivCoef( par = par, data = data3 ) ) )
      }

      # perform fit
      result$nls.lm <- nls.lm( par = startVal, fn = residFun, data = estData,
         jac = jac, ... )
      result$coefficients <- result$nls.lm$par
      result$iter <- result$nls.lm$niter
      result$convergence <- result$nls.lm$info > 0 && result$nls.lm$info < 5
      result$message <- result$nls.lm$message
   } else if( method == "Newton" ) {
      cesRss2 <- function( par, data ) {
         result <- cesRss( par = par, data = data )
         attributes( result )$gradient <- cesRssDeriv( par = par, data = data )
         return( result )
      }
      # save current setting for warning messages and suppress warning messages
      warnSaved <- options()$warn
      options( warn = -1 )
      # perform fit
      result$nlm <- nlm( f = cesRss2, p = startVal, data = estData, ... )
      # restore previous setting for warning messages
      options( warn = warnSaved )
      # extract results
      result$coefficients <- result$nlm$estimate
      result$iter <- result$nlm$iterations
      result$convergence <- result$nlm$code <= 2
   } else {
      stop( "argument 'method' must be either 'Nelder-Mead', 'BFGS',",
         " 'CG', 'L-BFGS-B', 'SANN', 'LM', 'Newton', or 'Kmenta'" )
   }

   # return also the call
   result$call <- matchedCall

   # returned the method used for the estimation
   result$method <- method

   # fitted values
   result$fitted.values <- cesCalc( xNames = xNames, data = data,
      coef = result$coefficients )

   # residuals
   result$residuals <- estData$y - result$fitted.values

   # unscaled covariance matrix
   gradients <- cesDerivCoef( par = result$coefficients, data = estData )
   result$cov.unscaled <- try( chol2inv( chol( crossprod( gradients ) ) ),
      silent = TRUE )
   if( !is.matrix( result$cov.unscaled ) ) {
      result$cov.unscaled <- matrix( NA, nrow = length( result$coefficients ),
         ncol = length( result$coefficients ) )
   }
   rownames( result$cov.unscaled ) <- names( result$coefficients )
   colnames( result$cov.unscaled ) <- names( result$coefficients )

   # nonlinear least squares
#    result$nls <- nls(
#       y ~ gamma * ( delta * x1^rho + ( 1 - delta ) * x2^rho )^(1/rho),
#       data = estData, start = result$startVal, trace = TRUE,
#       algorithm = "port", lower = c( -Inf, 0.01 , -Inf ),
#       upper = c( Inf, 0.99, Inf ) )

   class( result ) <- c( "cesEst", class( result ) )
   return( result )
}

