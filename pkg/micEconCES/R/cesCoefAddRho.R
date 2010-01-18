cesCoefAddRho <- function( coef, vrs, rho ) {

   if( !is.null( rho ) ) {
      if( vrs ) {
         coef <- c( coef[ -length( coef ) ], rho, coef[ length( coef ) ] )
      } else {
         coef <- c( coef, rho )
      }
      if( !is.null( names( coef ) ) ) {
         names( coef )[ length( coef ) - vrs ] <- "rho"
      }
   }
   return( coef )
}
