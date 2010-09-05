cesCalc <- function( xNames, data, coef, nested = FALSE, rhoApprox = 5e-6 ) {

   # check number of exogenous variables
   nExog <- length( xNames )
   if( nExog < 2 ) {
      stop( "argument 'xNames' must include the names of at least 2 variables" )
   }

   # check 
   if( nested && nExog != 4 ) {
      stop( "the nested CES function is currently implemented only for",
         " 4 inputs" )
   }

   # check names of exogenous variables
   checkNames( xNames, names( data ) )

   # check for VRS
   if( nExog == 2 ) {
      vrs <- length( coef ) >= 4
   } else {
      vrs <- length( coef ) - nExog  >= 3
   }

   # check number of coefficients
   if( nExog == 2 && length( coef ) != 3 + vrs ) {
      stop( "a CES function with 2 exogenous variables and",
         ifelse( vrs, " variable", " constant" ), " returns to scale",
         " must have ", 3 + vrs, " coefficients",
         " but you provided ", length( coef ), " coefficients" )
   } else if( nExog > 2 && !nested && length( coef ) != nExog + 2 + vrs ) {
      stop( "a non-nested CES function with ", nExog, " exogenous variables and",
         ifelse( vrs, " variable", " constant" ), " returns to scale",
         " must have ", nExog + 2 + vrs, " coefficients",
         " but you provided ", length( coef ), " coefficients" )
   } else if( nExog == 4 && nested && length( coef ) != nExog + 2 + vrs ) {
      stop( "a nested CES function with ", nExog, " exogenous variables and",
         ifelse( vrs, " variable", " constant" ), " returns to scale",
         " must have ", nExog + 2 + vrs, " coefficients",
         " but you provided ", length( coef ), " coefficients" )
   }

   # check for NAs in coefficients
   if( sum( is.na( coef ) ) > 0 ) {
      warning( "some of the coefficiencients are 'NA'" )
   }

   # names of coefficients
   coefNames <- cesCoefNames( nExog = nExog, vrs = vrs, nested = nested )

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
   if( !nested ) {
      deltaCoefs <- coef[ grep( "delta\\_", names( coef ) ) ]
      if( sum( is.na( deltaCoefs ) ) == 0 &&
            abs( sum( deltaCoefs, na.rm = TRUE ) - 1 ) / sum( abs( deltaCoefs ) ) >
            .Machine$double.eps^0.5 ) {
         stop( "the sum of the delta coefficients must sum up to 1" )
      }
   }

   # make the case of constant returns to scale (CRS) compatible to the VRS case
   if( !vrs ) {
      coef <- c( coef, nu = 1 )
   }

   # calculate the endogenous variable
   if( !nested ) {
      if( abs( coef[ "rho" ] ) <= rhoApprox ) {
         result <- log( coef[ "gamma" ] )
         for( i in 1:nExog ) {
            result <- result + coef[ paste( "delta", i, sep = "_" ) ] *
               coef[ "nu" ] * log( data[[ xNames[ i ] ]] )
         }
         for( i in 1:( nExog - 1 ) ) {
            for( j in ( i + 1 ):nExog ) {
               result <- result - 0.5 * coef[ "rho" ] * coef[ "nu" ] *
                  coef[ paste( "delta", i, sep = "_" ) ] *
                  coef[ paste( "delta", j, sep = "_" ) ] *
                  ( log( data[[ xNames[ i ] ]] ) - log( data[[ xNames[ j ] ]] ) )^2
            }
         }
         result <- exp( result )
      } else {
         if( coef[ "rho" ] == 0 ) {
            result <- NaN
         } else {
            result <- 0
            for( i in 1:nExog ) {
               result <- result + coef[ paste( "delta", i, sep = "_" ) ] *
                  data[[ xNames[ i ] ]]^( -coef[ "rho" ] )
            }
            result <- result^( -coef[ "nu" ] / coef[ "rho" ] )
            result <- coef[ "gamma" ] * result
         }
      }
   } else {
      result <- coef[ "gamma" ] * (
            ( coef[ "delta_1" ] * data[[ xNames[ 1 ] ]]^( -coef[ "rho_1" ] ) +
               ( 1 - coef[ "delta_1" ] ) * data[[ xNames[ 2 ] ]]^( -coef[ "rho_1" ] ) 
            )^( coef[ "rho" ] / coef[ "rho_1" ] ) +
            ( coef[ "delta_2" ] * data[[ xNames[ 3 ] ]]^( -coef[ "rho_2" ] ) +
               ( 1 - coef[ "delta_2" ] ) * data[[ xNames[ 4 ] ]]^( -coef[ "rho_2" ] ) 
            )^( coef[ "rho" ] / coef[ "rho_2" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] )
   }

   return( result )
}
