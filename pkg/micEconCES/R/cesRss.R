cesRss <- function( par, yName, xNames, data, vrs, rho = NULL,
      rhoApprox, nested = FALSE ) {

   # check rhoApprox
   if( !is.vector( rhoApprox ) || ! length( rhoApprox ) %in% c( 1, 5 ) ||
         !is.numeric( rhoApprox ) ) {
      stop( "argument 'rhoApprox' must be a numeric scalar",
         " or a numeric vector with exactly 5 elements" )
   }
   if( length( rhoApprox ) == 5 ) {
      rhoApprox <- rhoApprox[ 1 ]
   }

   # add coefficient 'rho' if it is fixed
   par <- cesCoefAddRho( coef = par, vrs = vrs, rho = rho )

   yHat <- cesCalc( xNames = xNames, data = data, coef = par,
      rhoApprox = rhoApprox, nested = nested )

   result <- sum( ( data[[ yName ]] - yHat )^2 )
   if( is.na( result ) ) {
      result <- Inf
   }
   return( result )
}
