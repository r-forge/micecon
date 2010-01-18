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
      phi <- par[ "phi" ]
   } else {
      phi <- 1
   }

   # derivatives with respect to gamma
   result[ , "gamma" ] <-
      ( delta * data[[ xNames[ 1 ] ]]^(-rho) + ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) )^( -phi / rho )

   # derivatives with respect to delta
   result[ , "delta" ] <- gamma * ( -phi / rho ) *
      ( delta * data[[ xNames[ 1 ] ]]^(-rho) + ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) )^( -phi / rho - 1 ) *
      ( data[[ xNames[ 1 ] ]]^(-rho) - data[[ xNames[ 2 ] ]]^(-rho) )

   # derivatives with respect to rho
   if( returnRho ) {
      result[ , "rho" ] <- gamma *
         log( delta * data[[ xNames[ 1 ] ]]^(-rho) + ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) ) *
         ( delta * data[[ xNames[ 1 ] ]]^(-rho) + ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) )^( -phi / rho ) *
         ( phi / rho^2 ) +
         gamma * ( phi / rho ) *
         ( delta * data[[ xNames[ 1 ] ]]^(-rho) + ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) )^( -phi / rho - 1 ) *
         ( delta * log( data[[ xNames[ 1 ] ]] ) * data[[ xNames[ 1 ] ]]^(-rho) +
            ( 1 - delta ) * log( data[[ xNames[ 2 ] ]] ) * data[[ xNames[ 2 ] ]]^(-rho) )
   }

   # derivatives with respect to phi
   if( vrs ) {
      result[ , "phi" ] <- gamma *
         log( delta * data[[ xNames[ 1 ] ]]^(-rho) + ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) ) *
         ( delta * data[[ xNames[ 1 ] ]]^(-rho) + ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) )^( -phi / rho ) *
         ( -1 / rho )
   }

   return( result )
}
