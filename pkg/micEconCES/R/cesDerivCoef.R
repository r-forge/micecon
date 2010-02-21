cesDerivCoef <- function( par, xNames, data, vrs, returnRho = TRUE,
      rhoApprox = 5e-6 ) {

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
   if( abs( rho ) > rhoApprox ) {
      result[ , "gamma" ] <-
         ( delta * data[[ xNames[ 1 ] ]]^(-rho) + ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) )^( -nu / rho )
   } else {
      result[ , "gamma" ] <- 
         data[[ xNames[ 1 ] ]]^( nu * delta ) *
         data[[ xNames[ 2 ] ]]^( nu * ( 1 - delta ) ) *
         exp( - 0.5 * rho * nu * delta * ( 1 - delta ) * 
         ( log( data[[ xNames[ 1 ] ]] ) - log( data[[ xNames[ 2 ] ]] ) )^2 )
   }

   # derivatives with respect to delta
   if( abs( rho ) > rhoApprox ) {
      result[ , "delta" ] <- - ( gamma * nu / rho ) *
         ( data[[ xNames[ 1 ] ]]^(-rho) - data[[ xNames[ 2 ] ]]^(-rho) ) *
         ( delta * data[[ xNames[ 1 ] ]]^(-rho) +
            ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) )^( - nu / rho - 1 )
   } else {
      result[ , "delta" ] <- gamma * nu *
         ( log( data[[ xNames[ 1 ] ]] ) - log( data[[ xNames[ 2 ] ]] ) ) *
         data[[ xNames[ 1 ] ]]^( nu * delta ) *
         data[[ xNames[ 2 ] ]]^( nu * ( 1 - delta ) ) *
         ( 1 - ( rho / 2 ) * ( 1 - 2 * delta + nu * delta * ( 1 - delta ) *
         ( log( data[[ xNames[ 1 ] ]] ) - log( data[[ xNames[ 2 ] ]] ) ) ) *
         ( log( data[[ xNames[ 1 ] ]] ) - log( data[[ xNames[ 2 ] ]] ) ) )
   }

   # derivatives with respect to rho
   if( returnRho ) {
      if( abs( rho ) > rhoApprox * 200 ) {
         result[ , "rho" ] <- ( gamma * nu / rho^2 ) *
            log( delta * data[[ xNames[ 1 ] ]]^(-rho) +
               ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) ) *
            ( delta * data[[ xNames[ 1 ] ]]^(-rho) +
               ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) )^( -nu / rho ) +
            ( gamma * nu / rho ) *
            ( delta * log( data[[ xNames[ 1 ] ]] ) * data[[ xNames[ 1 ] ]]^(-rho) +
               ( 1 - delta ) * log( data[[ xNames[ 2 ] ]] ) * data[[ xNames[ 2 ] ]]^(-rho) ) *
            ( delta * data[[ xNames[ 1 ] ]]^(-rho) +
               ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) )^( -nu / rho - 1 )
      } else {
         result[ , "rho" ] <- gamma * nu * delta * ( 1 - delta ) *
            data[[ xNames[ 1 ] ]]^( nu * delta ) *
            data[[ xNames[ 2 ] ]]^( nu * ( 1 - delta ) ) *
            ( - ( 1 / 2 ) * 
            ( log( data[[ xNames[ 1 ] ]] ) - log( data[[ xNames[ 2 ] ]] ) )^2
            + ( 1 / 3 ) * rho * ( 1 - 2 * delta ) *
            ( log( data[[ xNames[ 1 ] ]] ) - log( data[[ xNames[ 2 ] ]] ) )^3
            + ( 1 / 4 ) * rho * nu * delta * ( 1 - delta ) *
            ( log( data[[ xNames[ 1 ] ]] ) - log( data[[ xNames[ 2 ] ]] ) )^4 )
      }
   }

   # derivatives with respect to nu
   if( vrs ) {
      if( abs( rho ) > rhoApprox ) {
         result[ , "nu" ] <- - ( gamma / rho ) *
            log( delta * data[[ xNames[ 1 ] ]]^(-rho) +
               ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) ) *
            ( delta * data[[ xNames[ 1 ] ]]^(-rho) +
               ( 1 - delta ) * data[[ xNames[ 2 ] ]]^(-rho) )^( -nu / rho )
      } else {
         result[ , "nu" ] <- gamma * 
            data[[ xNames[ 1 ] ]]^( nu * delta ) *
            data[[ xNames[ 2 ] ]]^( nu * ( 1 - delta ) ) *
            ( delta * log( data[[ xNames[ 1 ] ]] ) + 
            ( 1 - delta ) * log( data[[ xNames[ 2 ] ]] ) -
            ( rho * delta * ( 1 - delta ) / 2 ) * 
            ( log( data[[ xNames[ 1 ] ]] ) - log( data[[ xNames[ 2 ] ]] ) )^2 *
            ( 1 + nu * ( delta * log( data[[ xNames[ 1 ] ]] ) + 
            ( 1 - delta ) * log( data[[ xNames[ 2 ] ]] ) ) ) ) 
      }
   }

   return( result )
}
