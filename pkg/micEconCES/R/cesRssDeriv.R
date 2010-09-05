cesRssDeriv <- function( par, yName, xNames, data, vrs, rho = NULL,
      rhoApprox = c( 5e-6, 5e-6, 5e-6, 1e-3, 5e-6 ), nested = FALSE ) {

   # check rhoApprox
   if( !is.vector( rhoApprox ) || length( rhoApprox ) != 5 ||
         !is.numeric( rhoApprox ) ) {
      stop( "argument 'rhoApprox' must be a numeric vector with exactly",
         "5 elements" )
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
