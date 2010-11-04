cesEstStart <- function( yName, xNames, data, vrs,
      method, start, rho, nParam, nested = FALSE ) {

   # number of explanatory variables
   nExog <- length( xNames )

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
         if( nested && nExog == 3 ) {
            start <- c( 1, 1, 0.5, 0.5, 0.25, rhoStart )
         } else if( nested && nExog == 4 ) {
            start <- c( 1, 0.5, 0.5, 0.5, 0.25, 0.25, rhoStart )
         } else if( !nested && nExog == 2 ) {
            start <- c( 1, 0.5, rhoStart )
         } else {
            stop( "cannot create starting values for a",
               ifelse( nested, " nested", " non-nested" ), " CES function",
               " with ", nExog, " explanatory variables" )
         }
         if( vrs ) {
            start <- c( start, 1 )
         }
         yTemp <- cesCalc( xNames = xNames, data = data, coef = start,
            nested = nested )
         start[ 1 + ( nested && nExog == 3 ) ] <- 
            mean( data[[ yName ]], na.rm = TRUE ) /
            mean( yTemp, na.rm = TRUE )
         if( !is.null( rho ) ) {
            start <- start[ -ifelse( nested, 3 + nExog, 3 ) ]
         }
      }
      if( length( start ) != nParam ) {
         stop( "wrong number of starting values:",
            " you provided ", length( start ), " values",
            " but the model has ", nParam, " parameters" )
      }
      names( start ) <- cesCoefNames( nExog = length( xNames ), vrs = vrs,
         returnRho = is.null( rho ), nested = nested )
      # checking starting values
      if( any( is.infinite( start ) ) ) {
         stop( "all starting values must be finite" )
      }
      # checking gamma
      if( nested && nExog == 3 ) {
         if( start[ "gamma_1" ] <= 0 ) {
            stop( "the starting value for 'gamma_1' must be positive" )
         } else if( start[ "gamma_2" ] <= 0 ) {
            stop( "the starting value for 'gamma_2' must be positive" )
         }
      } else {
         if( start[ "gamma" ] <= 0 ) {
            stop( "the starting value for 'gamma' must be positive" )
         }
      }
      # checking delta
      if( nested ) {
         if( start[ "delta_1" ] < 0 || start[ "delta_1" ] > 1 ) {
            stop( "the starting value for 'delta_1' must be between 0 and 1" )
         }
         if( start[ "delta_2" ] < 0 || start[ "delta_2" ] > 1 ) {
            stop( "the starting value for 'delta_2' must be between 0 and 1" )
         }
      } else {
         if( start[ "delta" ] < 0 || start[ "delta" ] > 1 ) {
            stop( "the starting value for 'delta' must be between 0 and 1" )
         }
      }
      # checking rho
      if( is.null( rho ) ) {
         if( start[ "rho" ] < -1 ) {
            stop( "the starting value for 'rho' must be -1 or larger" )
         }
      }
      if( nested ) {
         if( start[ "rho_1" ] < -1 ) {
            stop( "the starting value for 'rho_1' must be -1 or larger" )
         }
         if( nExog == 4 ) {
            if( start[ "rho_2" ] < -1 ) {
               stop( "the starting value for 'rho_2' must be -1 or larger" )
            }
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
