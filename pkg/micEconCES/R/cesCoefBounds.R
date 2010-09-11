cesCoefBounds <- function( vrs, returnRho, method, lower, nested = FALSE ) {

   if( method %in% c( "L-BFGS-B", "PORT", "DE" ) ) {
      if( lower ) {
         if( nested ) {
            result <- c( 0, 0, 0, -1, -1 )
         } else {
            result <- c( 0, 0 )
         }
      } else {
         if( nested ) {
            result <- c( Inf, 1, 1, Inf, Inf )
         } else {
            result <- c( Inf, 1 )
         }
      }
      if( returnRho ) {
         result <- c( result, ifelse( lower, -1, Inf ) )
      }
      if( vrs ) {
         result <- c( result, ifelse( lower, 0, Inf ) )
      }
   } else {
      result <- ifelse( lower, -Inf, Inf )
   }

   if( method == "DE" ) {
      result[ 1:2 ][ !is.finite( result[ 1:2 ] ) ] <- 1e10
      result[ !is.finite( result ) ] <- 10
   }

   return( result )
}
