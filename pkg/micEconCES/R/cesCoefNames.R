cesCoefNames <- function( nExog, vrs ) {

   if( nExog == 2 ) {
      result <- c( "gamma", "delta", "rho" )
   } else {
      result <- c( "gamma", paste( "delta", 1:nExog, sep = "_" ), "rho" )
   }
   if( vrs ) {
      result <- c( result, "phi" )
   }
   return( result )
}