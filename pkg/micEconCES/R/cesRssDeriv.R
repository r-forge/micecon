cesRssDeriv <- function( par, yName, xNames, data, vrs, rho = NULL,
      rhoApprox, nested = FALSE ) {

   # number of exogenous variables
   nExog <- length( xNames )

   # obtain names of coefficients
   coefNames <- cesCoefNames( nExog = nExog, vrs = vrs, 
      returnRho = is.null( rho ), nested = nested )

   # check rhoApprox
   if( !nested ) {
      cesCheckRhoApprox( rhoApprox = rhoApprox, elemNames = c( "y", coefNames ) )
   }

   # add coefficient 'rho' if it is fixed
   par <- cesCoefAddRho( coef = par, vrs = vrs, rho = rho )

   # calculate fitted values and residuals
   yHat <- cesCalc( xNames = xNames, data = data, coef = par,
      rhoApprox = rhoApprox[1], nested = nested )
   resid <- data[[ yName ]] - yHat

   # obtain derivatives of the CES with respect to coefficients
   derivCoef <- cesDerivCoef( par = par, xNames = xNames, data = data, 
      vrs = vrs, returnRho = is.null( rho ), rhoApprox = rhoApprox[-1],
      nested = nested )

   # prepare vector of gradients (to be returned)
   result <- numeric( ncol( derivCoef ) )
   names( result ) <- colnames( derivCoef )
   for( coefName in colnames( derivCoef ) ) {
      result[ coefName ] <- sum( - 2 * resid * derivCoef[ , coefName ] )
   }
   return( result )
}
