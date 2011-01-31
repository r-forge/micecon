cesCoefNames <- function( nExog, vrs, returnRho1 = TRUE, returnRho2 = TRUE, 
      returnRho = TRUE, nested = FALSE ) {

   if( nExog == 2 ) {
      result <- c( "gamma", "delta" )
   } else if( !nested ) {
      result <- c( "gamma", paste( "delta", 1:nExog, sep = "_" ) )
   } else if( nested && nExog == 3 ) {
      result <- c( "gamma", "delta_1", "delta_2" )
   } else if( nested && nExog == 4 ) {
      result <- c( "gamma", "delta_1", "delta_2", "delta_3" )
   } else {
      stop( "internal error: non-supported arguments to cesCoefNames()" )
   }
   if( returnRho1 && nested ) {
      result <- c( result, "rho_1" )
   }
   if( returnRho2 && nested && nExog == 4 ) {
      result <- c( result, "rho_2" )
   }
   if( returnRho ) {
      result <- c( result, "rho" )
   }
   if( vrs ) {
      result <- c( result, "nu" )
   }
   return( result )
}
