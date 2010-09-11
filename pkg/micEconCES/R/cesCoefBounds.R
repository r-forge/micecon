cesCoefBounds <- function( vrs, returnRho, method, lower, nExog, 
      nested = FALSE ) {

   if( method %in% c( "L-BFGS-B", "PORT", "DE" ) ) {
      if( lower ) {
         if( nested && nExog == 4 ) {
            result <- c( 0, 0, 0, -1, -1 )
         } else if( !nested && nExog == 2 ) {
            result <- c( 0, 0 )
         } else {
            stop( "cannot create lower bounds for the parameters of a",
               ifelse( nested, " nested", " non-nested" ), " CES function",
               " with ", nExog, " explanatory variables" )
         }
      } else {
         if( nested && nExog == 4 ) {
            result <- c( Inf, 1, 1, Inf, Inf )
         } else if( !nested && nExog == 2 ) {
            result <- c( Inf, 1 )
         } else {
            stop( "cannot create upper bounds for the parameters of a",
               ifelse( nested, " nested", " non-nested" ), " CES function",
               " with ", nExog, " explanatory variables" )
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
