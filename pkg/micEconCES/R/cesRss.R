cesRss <- function( par, yName, xNames, data, vrs, rho = NULL,
      rhoApprox = 5e-6 ) {

   # add coefficient 'rho' if it is fixed
   par <- cesCoefAddRho( coef = par, vrs = vrs, rho = rho )

   yHat <- cesCalc( xNames = xNames, data = data, coef = par,
      rhoApprox = rhoApprox )

   result <- sum( ( data[[ yName ]] - yHat )^2 )
   if( is.na( result ) ) {
      result <- Inf
   }
   return( result )
}
