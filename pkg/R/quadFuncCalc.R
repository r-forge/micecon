quadFuncCalc <- function( xNames, data, coef, shifterNames = NULL,
      quadHalf = TRUE ) {

   checkNames( c( xNames, shifterNames ), names( data ) )

   nExog <- length( xNames )
   nShifter <- length( shifterNames )
   nCoef <- 1 + nExog + nExog * ( nExog + 1 ) / 2 + nShifter

   if( nCoef != length( coef ) ) {
      stop( "a quadratic function with ", nExog, " exogenous variables",
         " and ", nShifter, "shifter variables",
         " must have exactly ", nCoef, " coefficients" )
   }

   result <- rep( coef[ "a_0" ], nrow( data ) )
   for( i in seq( along = xNames ) ) {
      result <- result + coef[ paste( "a", i, sep = "_" ) ] * 
         data[[ xNames[ i ] ]]
      for( j in seq( along = xNames ) ) {
         result <- result + ifelse( quadHalf, 0.5, 1 ) * 
            coef[ paste( "b", min( i, j ), max( i, j ), sep = "_" ) ] *
            data[[ xNames[ i ] ]]  * data[[ xNames[ j ] ]]
      }
   }
   for( i in seq( along = shifterNames ) ) {
      result <- result + coef[ paste( "d", i, sep = "_" ) ] * 
         data[[ shifterNames[ i ] ]]
   }

   names( result ) <- rownames( data )
   return( result )
}
