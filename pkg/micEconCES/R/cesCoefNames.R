cesCoefNames <- function( nExog, vrs ) {
   if( nExog != 2 ) {
      stop( "currently coefficient names can be created",
         " only for CES functions with exactly 2 exogenous variables" )
   }
   result <- c( "gamma", "delta", "rho" )
   if( vrs ) {
      result <- c( result, "phi" )
   }
   return( result )
}