cesRssDeriv <- function( par, data ) {
   result <- par

   yHat <- cesCalc( xNames = names( data )[ -1 ], data = data, coef = par )
   resid <- data$y - yHat

   derivCoef <- cesDerivCoef( par = par, data = data )

   for( coefName in colnames( derivCoef ) ) {
      result[ coefName ] <- sum( - 2 * resid * derivCoef[ , coefName ] )
   }

   return( result )
}
