quadFuncCalc <- function( xNames, data, coef, quadHalf = TRUE ) {

   checkNames( c( xNames ), names( data ) )

   nExog <- length( xNames )
   nCoef <- 1 + nExog + nExog * ( nExog + 1 ) / 2

   if( nCoef != length( coef ) ) {
      stop( "a quadratic function with ", nExog, " exogenous variables",
         " must have exactly ", nCoef, " coefficients" )
   }

   alpha0 <- coef[ "a_0" ]
   alpha <- rep( NA, nExog )
   for( i in seq( along = xNames ) ) {
      alpha[ i ] <- coef[ paste( "a", i, sep = "_" ) ]
   }
   beta <- matrix( NA, nrow = nExog, ncol = nExog )
   for( i in seq( along = xNames ) ) {
      for( j in seq( along = xNames ) ) {
         beta[ i, j ] <- coef[ paste( "b", min( i, j ), max( i, j ), sep = "_" ) ]
      }
   }

   result <- rep( alpha0, nrow( data ) )
   for( i in 1:nExog ) {
      result <- result + alpha[ i ] * data[[ xNames[ i ] ]]
      for( j in 1:nExog ) {
         result <- result + ifelse( quadHalf, 0.5, 1 ) * beta[ i, j ] *
            data[[ xNames[ i ] ]]  * data[[ xNames[ j ] ]]
      }
   }

   names( result ) <- rownames( data )
   return( result )
}
