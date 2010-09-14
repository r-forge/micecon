cesDerivCoef <- function( par, xNames, data, vrs, nested = FALSE, 
      returnRho = TRUE, rhoApprox = c( 5e-6, 5e-6, 1e-3, 5e-6 ) ) {

   # check rhoApprox
   if( !is.vector( rhoApprox ) || length( rhoApprox ) != 4 ||
         !is.numeric( rhoApprox ) ) {
      stop( "argument 'rhoApprox' must be a numeric vector with exactly",
         "4 elements" )
   }

   # number of exogenous variables
   nExog <- length( xNames )

   # names of coefficients
   coefNames <- cesCoefNames( nExog = nExog, vrs = vrs, returnRho = returnRho,
      nested = nested )

   # derivatives of the CES with respect to the coefficients/parameters
   result <- matrix( NA, nrow = nrow( data ), ncol = length( coefNames ) )
   colnames( result ) <- coefNames
   names( par ) <- cesCoefNames( nExog = nExog, vrs = vrs, returnRho = TRUE,
      nested = nested )

   if( !nested ) {
      if( nExog != 2 ) {
         stop( "the derivatives of the non-nested CES can be calculated",
            " only for two inputs" )
      }
      gamma <- par[ "gamma" ]
      delta <- par[ "delta" ]
      rho <- par[ "rho" ]
      if( vrs ) {
         nu <- par[ "nu" ]
      } else {
         nu <- 1
      }

      # derivatives with respect to gamma
      if( abs( rho ) > rhoApprox[1] ) {
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
      if( abs( rho ) > rhoApprox[2] ) {
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
         if( abs( rho ) > rhoApprox[3] ) {
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
         if( abs( rho ) > rhoApprox[4] ) {
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
   } else if( nExog == 3 ) { # nested CES with 3 outputs
      gamma1 <- par[ "gamma_1" ]
      gamma2 <- par[ "gamma_2" ]
      delta1 <- par[ "delta_1" ]
      delta2 <- par[ "delta_2" ]
      rho1 <- par[ "rho_1" ]
      rho <- par[ "rho" ]
      if( vrs ) {
         nu <- par[ "nu" ]
      } else {
         nu <- 1
      }

      # main parts of the nested CES with 3 inputs
      B1 <- delta1 * data[[ xNames[ 1 ] ]]^(-rho1) + 
         ( 1 - delta1 ) * data[[ xNames[ 2 ] ]]^(-rho1)
      B <- delta2 * gamma1^(-rho) * B1^( rho / rho1 ) + 
         ( 1 - delta2 ) * data[[ xNames[ 3 ] ]]^(-rho)
    
      # derivatives with respect to gamma_1 and gamma_2
      result[ , "gamma_1" ] <- gamma2 * nu * delta2 * gamma1^(-rho-1) * 
         B1^(rho/rho1) * B^(-nu/rho-1)
      result[ , "gamma_2" ] <- B^(-nu/rho)

      # derivatives with respect to delta_1 and delta_2
      result[ , "delta_1" ] <- - gamma2 * nu * delta2 * gamma1^(-rho) * 
         B1^(rho/rho1-1) / rho1 * B^(-nu/rho-1) *
         ( data[[ xNames[ 1 ] ]]^(-rho1) - data[[ xNames[ 2 ] ]]^(-rho1) )
      result[ , "delta_2" ] <- - gamma2 * ( nu / rho ) *
         ( gamma1^(-rho) * B1^(rho/rho1) - data[[ xNames[ 3 ] ]]^(-rho) ) *
         B^(-nu/rho-1)

      # derivatives with respect to rho_1
      result[ , "rho_1" ] <- - gamma2 * nu * delta2 / rho1 * gamma1^(-rho) * 
         B1^(rho/rho1) * B^(-nu/rho-1) *
         ( - log( B1 ) / rho1 +
         ( - delta1 * log( data[[ xNames[ 1 ] ]] ) * data[[ xNames[ 1 ] ]]^(-rho1) -
         ( 1 - delta1 ) * log( data[[ xNames[ 2 ] ]] ) *
            data[[ xNames[ 2 ] ]]^(-rho1) ) / B1 )

      # derivatives with respect to rho
      if( returnRho ) {
         result[ , "rho" ] <- - gamma2  * ( nu / rho ) * B^(-nu/rho) *
            ( -log( B ) /rho +
               B^(-1) * ( delta2 * gamma1^(-rho) * B1^(rho/rho1) *
                  ( - log( gamma1 ) + log( B1 ) / rho1 ) -
               ( 1 - delta2 ) * log( data[[ xNames[ 3 ] ]] ) *
               data[[ xNames[ 3 ] ]]^(-rho) ) )
      }

      # derivatives with respect to nu
      if( vrs ) {
         result[ , "nu" ] <- - gamma2 * log( B ) * B^( -nu / rho ) / rho
      }

   } else if( nExog == 4 ) { # nested CES with 4 outputs
      gamma <- par[ "gamma" ]
      delta1 <- par[ "delta_1" ]
      delta2 <- par[ "delta_2" ]
      rho1 <- par[ "rho_1" ]
      rho2 <- par[ "rho_2" ]
      rho <- par[ "rho" ]
      if( vrs ) {
         nu <- par[ "nu" ]
      } else {
         nu <- 1
      }

      # main parts of the nested CES with 4 inputs
      B1 <- delta1 * data[[ xNames[ 1 ] ]]^(-rho1) + 
         ( 1 - delta1 ) * data[[ xNames[ 2 ] ]]^(-rho1)
      B2 <- delta2 * data[[ xNames[ 3 ] ]]^(-rho2) + 
         ( 1 - delta2 ) * data[[ xNames[ 4 ] ]]^(-rho2)
      B <- B1^( rho / rho1 ) + B2^( rho / rho2 )
    
      # derivatives with respect to gamma
      result[ , "gamma" ] <- B^(-nu/rho)

      # derivatives with respect to delta_1 and delta_2
      result[ , "delta_1" ] <- gamma * ( -nu / rho ) * B^((-nu-rho)/rho) *
         (rho/rho1) * B1^((rho-rho1)/rho1) * 
         ( data[[ xNames[ 1 ] ]]^(-rho1) - data[[ xNames[ 2 ] ]]^(-rho1) )
      result[ , "delta_2" ] <- gamma * ( -nu / rho ) * B^((-nu-rho)/rho) *
         (rho/rho2) * B2^((rho-rho2)/rho2) * 
         ( data[[ xNames[ 3 ] ]]^(-rho2) - data[[ xNames[ 4 ] ]]^(-rho2) )

      # derivatives with respect to rho_1 and rho_2
      result[ , "rho_1" ] <- gamma * ( -nu / rho ) * B^((-nu-rho)/rho) *
         ( log( B1 ) * B1^(rho/rho1) * ( -rho/rho1^2 ) + 
         B1^((rho-rho1)/rho1) * (rho/rho1) *
         ( -delta1 * log( data[[ xNames[ 1 ] ]] ) * data[[ xNames[ 1 ] ]]^(-rho1) - 
         ( 1 - delta1 ) * log( data[[ xNames[ 2 ] ]] ) * data[[ xNames[ 2 ] ]]^(-rho1) ) )
      result[ , "rho_2" ] <- gamma * ( -nu / rho ) * B^((-nu-rho)/rho) *
         ( log( B2 ) * B2^(rho/rho2) * ( -rho/rho2^2 ) + 
         B2^((rho-rho2)/rho2) * (rho/rho2) *
         ( -delta2 * log( data[[ xNames[ 3 ] ]] ) * data[[ xNames[ 3 ] ]]^(-rho2) - 
         ( 1 - delta2 ) * log( data[[ xNames[ 4 ] ]] ) * data[[ xNames[ 4 ] ]]^(-rho2) ) )

      # derivatives with respect to rho
      if( returnRho ) {
         result[ , "rho" ] <- gamma * log( B ) * B^(-nu/rho) * ( nu / rho^2 ) +
            gamma * ( -nu/rho ) * B^((-nu-rho)/rho) * 
            ( log( B1 ) * B1^(rho/rho1) / rho1 + log( B2 ) * B2^(rho/rho2) / rho2 )
      }

      # derivatives with respect to nu
      if( vrs ) {
         result[ , "nu" ] <- - gamma * log( B ) * B^( -nu / rho ) / rho
      }
   } else {
      stop( "the derivatives of the nested CES can be calculated",
         " only for three and four inputs" )
   }

   return( result )
}
