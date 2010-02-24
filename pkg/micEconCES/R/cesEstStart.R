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
      # checking starting values
      if( any( is.infinite( start ) ) ) {
         stop( "all starting values must be finite" )
      }
      # checking gamma
      if( start[ "gamma" ] <= 0 ) {
         stop( "the starting value for 'gamma' must be positive" )
      }
      # checking delta
      if( start[ "delta" ] < 0 || start[ "delta" ] > 1 ) {
         stop( "the starting value for 'delta' must be between 0 and 1" )
      }
      # checking rho
      if( is.null( rho ) ) {
         if( start[ "rho" ] < -1 ) {
            stop( "the starting value for 'rho' must be -1 or larger" )
         }
      }
      # checking nu
      if( vrs ) {
         if( start[ "nu" ] <= 0 ) {
            stop( "the starting value for 'nu' must be positive" )
         }
      }
   }

   return( start )
}
