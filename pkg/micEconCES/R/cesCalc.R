cesCalc <- function( xNames, data, coef ) {

   # check number of exogenous variables
   nExog <- length( xNames )
   if( nExog < 2 ) {
      stop( "argument 'xNames' must include the names of at least 2 variables" )
   }

   # check names of exogenous variables
   checkNames( xNames, names( data ) )

   # check number of coefficients
   if( nExog == 2 && ( length( coef ) < 3 || length( coef ) > 4 ) ) {
      stop( "a CES function with 2 exogenous variables",
         " must have either 3 (CRS) or 4 (VRS) coefficients" )
   } else if( nExog > 2 &&
         ( length( coef ) < nExog + 2 | length( coef ) > nExog + 3 ) ) {
      stop( "a CES function with ", nExog, " exogenous variables",
         " must have either ", nExog + 2, " (CRS) or ",
         nExog + 3, " (VRS) coefficients" )
   }

   # names of coefficients
   if( nExog == 2 ) {
      coefNames <- c( "gamma", "delta", "rho", "phi" )[ 1:length( coef ) ]
   } else {
      coefNames <- c( "gamma", paste( "delta", 1:nExog, sep = "_" ),
         "rho", "phi" )[ 1:length( coef ) ]
   }

   # assign or check names of coefficients
   if( is.null( names( coef ) ) ) {
      names( coef ) <- coefNames
   } else {
      coefTest <- coefNames %in% names( coef )
      if( any( !coefTest ) ) {
         stop( "following coefficient name(s) is/are missing in argument",
            " 'coef': ", paste( coefNames[ !coefTest ], collapse = ", " ) )
      }
   }

   # make the case of two explanatory compatible to the case of N variables
   if( nExog == 2 ) {
      names( coef )[ names( coef ) == "delta" ] <- "delta_1"
      coef <- c( coef, delta_2 = 1 - unname( coef[ "delta_1" ] ) )
   }

   # check if the deltas sum up to one
   if( abs( sum( coef[ grep( "delta\\_", names( coef ) ) ] ) - 1 ) >
         .Machine$double.eps ) {
      stop( "the sum of the delta coefficients must sum up to 1" )
   }

   # make the case of constant returns to scale (CRS) compatible to the VRS case
   if( ! "phi" %in% names( coef ) ) {
      coef <- c( coef, phi = 1 )
   }

   # calculate the endogenous variable
   result <- 0
   for( i in 1:nExog ) {
      result <- result + coef[ paste( "delta", i, sep = "_" ) ] *
         data[[ xNames[ i ] ]]^( -coef[ "rho" ] )
   }
   result <- result^( -coef[ "phi" ] / coef[ "rho" ] )
   result <- coef[ "gamma" ] * result 

   return( result )
}
