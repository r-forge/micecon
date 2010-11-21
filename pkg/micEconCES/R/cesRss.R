cesRss <- function( par, yName, xNames, data, vrs, rho1 = NULL, 
      rho2 = NULL, rho = NULL, rhoApprox, nested = FALSE ) {

   # check rhoApprox
   
   if( !is.vector( rhoApprox ) || !is.numeric( rhoApprox ) ) {
      stop( "argument 'rhoApprox' must be a numeric scalar",
         " or a numeric vector" )
   }
   rhoApprox <- cesCheckRhoApprox( rhoApprox = rhoApprox, withY = TRUE,
      withDeriv = NA )
   rhoApprox <- rhoApprox[ "y" ]

   # add coefficients 'rho_1', 'rho_2', and 'rho' if they are fixed
   par <- cesCoefAddRho( coef = par, vrs = vrs, rho1 = rho1, rho2 = rho2, 
      rho = rho, nExog = length( xNames ), nested = nested )

   yHat <- cesCalc( xNames = xNames, data = data, coef = par,
      rhoApprox = rhoApprox, nested = nested )

   result <- sum( ( data[[ yName ]] - yHat )^2 )
   if( is.na( result ) ) {
      result <- Inf
   }
   return( result )
}
