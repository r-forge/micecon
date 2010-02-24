cesEstStart <- function( yName, xNames, data, vrs,
      method, start, rho, nParam ) {

   # start values
   if( method %in% c( "Kmenta", "DE" ) ) {
      if( !is.null( start ) ) {
         warning( "ignoring starting values because they are not required",
            " for method '", method, "'" )
         start <- NULL
      }
   } else {
      if( is.null( start ) ) {
         rhoStart <- ifelse( is.null( rho ), 0.25, rho )
         start <- c( 1, 0.5, rhoStart, 1 )[ 1:( 3 + vrs ) ]
         yTemp <- cesCalc( xNames = xNames, data = data, coef = start )
         start[ 1 ] <- mean( data[[ yName ]], na.rm = TRUE ) /
            mean( yTemp, na.rm = TRUE )
         if( !is.null( rho ) ) {
            start <- start[ -3 ]
         }
      }
      if( length( start ) != nParam ) {
         stop( "wrong number of starting values:",
            " you provided ", length( start ), " values",
            " but the model has ", nParam, " parameters" )
      }
      names( start ) <- cesCoefNames( nExog = length( xNames ), vrs = vrs,
         returnRho = is.null( rho ) )
   }

   return( start )
}
