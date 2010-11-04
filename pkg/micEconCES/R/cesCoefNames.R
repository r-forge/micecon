cesCoefNames <- function( nExog, vrs, returnRho = TRUE, nested = FALSE ) {

   if( nExog == 2 ) {
      result <- c( "gamma", "delta" )
   } else if( !nested ) {
      result <- c( "gamma", paste( "delta", 1:nExog, sep = "_" ) )
   } else if( nested && nExog == 3 ) {
      result <- c( "gamma_1", "gamma_2", "delta_1", "delta_2", "rho_1" )
   } else if( nested && nExog == 4 ) {
      result <- c( "gamma", "delta_1", "delta_2", "delta_3", "rho_1", "rho_2" )
   } else {
      stop( "internal error: non-supported arguments to cesCoefNames()" )
   }
   if( returnRho ) {
      result <- c( result, "rho" )
   }
   if( vrs ) {
      result <- c( result, "nu" )
   }
   return( result )
}
