cesCoefBounds <- function( vrs, returnRho, method, lower, nested = FALSE ) {

   if( method %in% c( "L-BFGS-B", "PORT" ) ) {
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
   } else if( method == "DE" ) {
      if( lower ) {
         if( nested ) {
            result <- c( 0, 0, 0, -1, -1 )
         } else {
            result <- c( 0, 0 )
         }
      } else {
         if( nested ) {
            result <- c( 1e10, 1, 1, 10, 10 )
         } else {
            result <- c( 1e10, 1 )
         }
      }
      if( returnRho ) {
         result <- c( result, ifelse( lower, -1, 10 ) )
      }
      if( vrs ) {
         result <- c( result, ifelse( lower, 0, 10 ) )
      }
   } else {
      result <- ifelse( lower, -Inf, Inf )
   }
   return( result )
}
