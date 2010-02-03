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
   if( rho != 0 ) {
      result[ , "gamma" ] <-
         ( delta * data[[ xNames[ 1 ] ]]^(-rho) + ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) )^( -nu / rho )
   } else {
      result[ , "gamma" ] <- 
         data[[ xNames[ 1 ] ]]^( nu * delta ) *
         data[[ xNames[ 2 ] ]]^( nu * ( 1 - delta ) )
   }

   # derivatives with respect to delta
   if( rho != 0 ) {
      result[ , "delta" ] <- gamma * ( -nu / rho ) *
         ( delta * data[[ xNames[ 1 ] ]]^(-rho) + ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) )^( -nu / rho - 1 ) *
         ( data[[ xNames[ 1 ] ]]^(-rho) - data[[ xNames[ 2 ] ]]^(-rho) )
   } else {
      result[ , "delta" ] <- gamma * nu *
         ( log( data[[ xNames[ 1 ] ]] ) - log( data[[ xNames[ 2 ] ]] ) ) *
         data[[ xNames[ 1 ] ]]^( nu * delta ) *
         data[[ xNames[ 2 ] ]]^( nu * ( 1 - delta ) )
   }

   # derivatives with respect to rho
   if( returnRho ) {
      if( rho != 0 ) {
         result[ , "rho" ] <- gamma *
            log( delta * data[[ xNames[ 1 ] ]]^(-rho) + ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) ) *
            ( delta * data[[ xNames[ 1 ] ]]^(-rho) + ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) )^( -nu / rho ) *
            ( nu / rho^2 ) +
            gamma * ( nu / rho ) *
            ( delta * data[[ xNames[ 1 ] ]]^(-rho) + ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) )^( -nu / rho - 1 ) *
            ( delta * log( data[[ xNames[ 1 ] ]] ) * data[[ xNames[ 1 ] ]]^(-rho) +
               ( 1 - delta ) * log( data[[ xNames[ 2 ] ]] ) * data[[ xNames[ 2 ] ]]^(-rho) )
      } else {
         result[ , "nu" ] <- gamma * 
            ( delta * log( data[[ xNames[ 1 ] ]] ) + 
            ( 1 - delta ) * log( data[[ xNames[ 2 ] ]] ) ) *
            data[[ xNames[ 1 ] ]]^( nu * delta ) *
            data[[ xNames[ 2 ] ]]^( nu * ( 1 - delta ) )
      }
   }

   # derivatives with respect to nu
   if( vrs ) {
      if( rho != 0 ) {
         result[ , "nu" ] <- gamma *
            log( delta * data[[ xNames[ 1 ] ]]^(-rho) + ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) ) *
            ( delta * data[[ xNames[ 1 ] ]]^(-rho) + ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) )^( -nu / rho ) *
            ( -1 / rho )
      } else {
         result[ , "rho" ] <- - 0.5 * gamma * nu * delta * ( 1 - delta ) *
            data[[ xNames[ 1 ] ]]^( nu * delta ) *
            data[[ xNames[ 2 ] ]]^( nu * ( 1 - delta ) ) *
            ( log( data[[ xNames[ 1 ] ]] ) - log( data[[ xNames[ 2 ] ]] ) )^2
      }
   }

   return( result )
}
