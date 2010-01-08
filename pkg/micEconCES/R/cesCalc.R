cesCalc <- function( xNames, data, coef ) {

   checkNames( xNames, names( data ) )

   nExog <- length( xNames )

   if( length( coef ) < nExog + 1 | length( coef ) > nExog + 2 ) {
      stop( "a CES function with ", nExog, " exogenous variables",
         " must have either ", nExog + 1, " (CRS) or ",
         nExog + 2, " (VRS) coefficients" )
   }

   coefNames <- c( "gamma", "alpha", "rho", "phi" )[ 1:length( coef ) ]
   if( is.null( names( coef ) ) ) {
      names( coef ) <- coefNames
   } else {
      coefTest <- coefNames %in% names( coef )
      if( any( !coefTest ) ) {
         stop( "following coefficient name(s) is/are missing in argument",
            " 'coef': ", paste( coefNames[ !coefTest ], collapse = ", " ) )
      }
   }

   if( ! "phi" %in% names( coef ) ) {
      coef <- c( coef, phi = 1 )
   }

   result <- coef[ "gamma" ] *
      ( coef[ "alpha" ] * data[[ xNames[ 1 ] ]]^coef[ "rho" ] +
      ( 1 - coef[ "alpha" ] ) * data[[ xNames[ 2 ] ]]^coef[ "rho" ] )^
      ( coef[ "phi" ] / coef[ "rho" ] )

   return( result )
}