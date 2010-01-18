cesDerivCoef <- function( par, xNames, data, vrs, returnRho = TRUE ) {

   # names of coefficients
   coefNames <- cesCoefNames( nExog = 2, vrs = vrs, returnRho = returnRho )

   # derivatives of the CES with respect to the coefficients/parameters
   result <- matrix( NA, nrow = nrow( data ), ncol = length( coefNames ) )
   colnames( result ) <- coefNames
   names( par ) <- cesCoefNames( nExog = 2, vrs = vrs, returnRho = TRUE )

   gamma <- par[ "gamma" ]
   delta <- par[ "delta" ]
   rho <- par[ "rho" ]
   if( vrs ) {
      nu <- par[ "nu" ]
   } else {
      nu <- 1
   }

   # derivatives with respect to gamma
   result[ , "gamma" ] <-
      ( delta * data[[ xNames[ 1 ] ]]^(-rho) + ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) )^( -nu / rho )

   # derivatives with respect to delta
   result[ , "delta" ] <- gamma * ( -nu / rho ) *
      ( delta * data[[ xNames[ 1 ] ]]^(-rho) + ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) )^( -nu / rho - 1 ) *
      ( data[[ xNames[ 1 ] ]]^(-rho) - data[[ xNames[ 2 ] ]]^(-rho) )

   # derivatives with respect to rho
   if( returnRho ) {
      result[ , "rho" ] <- gamma *
         log( delta * data[[ xNames[ 1 ] ]]^(-rho) + ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) ) *
         ( delta * data[[ xNames[ 1 ] ]]^(-rho) + ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) )^( -nu / rho ) *
         ( nu / rho^2 ) +
         gamma * ( nu / rho ) *
         ( delta * data[[ xNames[ 1 ] ]]^(-rho) + ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) )^( -nu / rho - 1 ) *
         ( delta * log( data[[ xNames[ 1 ] ]] ) * data[[ xNames[ 1 ] ]]^(-rho) +
            ( 1 - delta ) * log( data[[ xNames[ 2 ] ]] ) * data[[ xNames[ 2 ] ]]^(-rho) )
   }

   # derivatives with respect to nu
   if( vrs ) {
      result[ , "nu" ] <- gamma *
         log( delta * data[[ xNames[ 1 ] ]]^(-rho) + ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) ) *
         ( delta * data[[ xNames[ 1 ] ]]^(-rho) + ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) )^( -nu / rho ) *
         ( -1 / rho )
   }

   return( result )
}
