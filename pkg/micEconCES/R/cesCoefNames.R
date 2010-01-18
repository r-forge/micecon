cesCoefNames <- function( nExog, vrs, returnRho = TRUE ) {

   if( nExog == 2 ) {
      result <- c( "gamma", "delta" )
   } else {
      result <- c( "gamma", paste( "delta", 1:nExog, sep = "_" ) )
   }
   if( returnRho ) {
      result <- c( result, "rho" )
   }
   if( vrs ) {
      result <- c( result, "phi" )
   }
   return( result )
}