quadFuncEla <- function( xNames, data, coef,
   yName = NULL, shifterNames = NULL, quadHalf = TRUE ) {

   checkNames( c( xNames, yName ), names( data ) )

   nExog <- length( xNames )
   nCoef <- 1 + nExog + nExog * ( nExog + 1 ) / 2

   if( nCoef > length( coef ) ) {
      stop( "a quadratic function with ", nExog, " exogenous variables",
         " must have at least ", nCoef, " coefficients" )
   }

   if( is.null( yName ) ){
      yHat <- quadFuncCalc( xNames = xNames, data = data, coef = coef, 
         shifterNames = shifterNames, quadHalf = quadHalf )
   } else {
      yHat <- data[[ yName ]]
   }

   result <- quadFuncDeriv( xNames = xNames, data = data, coef = coef,
      quadHalf = quadHalf )$deriv
   for( i in xNames ) {
      result[[ i ]] <- result[[ i ]] * data[[ i ]] / yHat
   }

   class( result ) <- c( "quadFuncEla", class( result ) )
   return( result )
}
