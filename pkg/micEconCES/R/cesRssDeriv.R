cesRssDeriv <- function( par, data ) {

   # calculate fitted values and residuals
   yHat <- cesCalc( xNames = names( data )[ -1 ], data = data, coef = par )
   resid <- data$y - yHat

   # obtain derivatives of the CES with respect to coefficients
   derivCoef <- cesDerivCoef( par = par, data = data )

   # prepare vector of gradients (to be returned)
   result <- numeric( length( par ) )
   names( result ) <- colnames( derivCoef )
   for( coefName in colnames( derivCoef ) ) {
      result[ coefName ] <- sum( - 2 * resid * derivCoef[ , coefName ] )
   }
   return( result )
}
