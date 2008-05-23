translogCalc <- function( xNames, data, allCoef, quadHalf = TRUE,
   dataLogged = FALSE ) {

   checkNames( c( xNames ), names( data ) )

   nExog <- length( xNames )
   nCoef <- 1 + nExog + nExog * ( nExog + 1 ) / 2

   if( nCoef != length( allCoef ) ) {
      stop( "a translog function with ", nExog, " exogenous variables",
         " must have exactly ", nCoef, " coefficients" )
   }

   if( dataLogged ) {
      logData <- data
   } else {
      logData <- data.frame( no = c( 1:nrow( data ) ) )
      for( i in seq( along = xNames ) ) {
         logData[[ xNames[ i ] ]] <- log( data[[ xNames[ i ] ]] )
      }
   }

   result <- quadFuncCalc( xNames, logData, allCoef, quadHalf = quadHalf )

   if( !dataLogged ) {
      result <- exp( result )
   }

   return( result )
}
