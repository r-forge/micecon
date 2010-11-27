cesDerivCoef <- function( par, xNames, data, vrs, nested = FALSE, 
      returnRho1 = TRUE, returnRho2 = TRUE, returnRho = TRUE, rhoApprox ) {

   # number of exogenous variables
   nExog <- length( xNames )

   # names of coefficients
   coefNames <- cesCoefNames( nExog = nExog, vrs = vrs, 
      returnRho1 = returnRho1, returnRho2 = returnRho2, 
      returnRho = returnRho, nested = nested )

   # check rhoApprox
   if( !nested ) {
      rhoApprox <- cesCheckRhoApprox( rhoApprox = rhoApprox, withY = NA,
         withDeriv = TRUE )
   }

   # derivatives of the CES with respect to the coefficients/parameters
   result <- matrix( NA, nrow = nrow( data ), ncol = length( coefNames ) )
   colnames( result ) <- coefNames
   names( par ) <- cesCoefNames( nExog = nExog, vrs = vrs, returnRho = TRUE,
      returnRho1 = TRUE, returnRho2 = TRUE, nested = nested )

   ###########################   non-nested CES   ##############################
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
      if( abs( rho ) > rhoApprox[ "gamma" ] ) {
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
      if( abs( rho ) > rhoApprox[ "delta" ] ) {
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
         if( abs( rho ) > rhoApprox[ "rho" ] ) {
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
         if( abs( rho ) > rhoApprox[ "nu" ] ) {
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

   ########################   nested CES with 3 inputs   #######################
   } else if( nExog == 3 ) { 
      if( !vrs ) {
         par <- c( par, nu = 1 )
      }

      # derivatives with respect to gamma_1
      result[ , "gamma_1" ] <- cesInterN3( 
         funcName = "cesDerivCoefN3Gamma1", par = par, 
         xNames = xNames, data = data, rhoApprox = rhoApprox[ "gamma" ] )

      # derivatives with respect to gamma_2
      result[ , "gamma_2" ] <- cesInterN3( 
         funcName = "cesDerivCoefN3Gamma2", par = par, 
         xNames = xNames, data = data, rhoApprox = rhoApprox[ "gamma" ] )

      # derivatives with respect to delta_1
      result[ , "delta_1" ] <- cesInterN3( 
         funcName = "cesDerivCoefN3Delta1", par = par, 
         xNames = xNames, data = data, rhoApprox = rhoApprox[ "delta" ] )

      # derivatives with respect to delta_2
      result[ , "delta_2" ] <- cesInterN3( 
         funcName = "cesDerivCoefN3Delta2", par = par, 
         xNames = xNames, data = data, rhoApprox = rhoApprox[ "delta" ] )

      # derivatives with respect to rho_1
      if( returnRho1 ) {
         result[ , "rho_1" ] <- cesInterN3( 
            funcName = "cesDerivCoefN3Rho1", par = par, 
            xNames = xNames, data = data, rhoApprox = rhoApprox[ "rho" ] )
      }

      # derivatives with respect to rho
      if( returnRho ) {
         result[ , "rho" ] <- cesInterN3( 
            funcName = "cesDerivCoefN3Rho", par = par, 
            xNames = xNames, data = data, rhoApprox = rhoApprox[ "rho" ] )
      }

      # derivatives with respect to nu
      if( vrs ) {
         result[ , "nu" ] <- cesInterN3( 
            funcName = "cesDerivCoefN3Nu", par = par, 
            xNames = xNames, data = data, rhoApprox = rhoApprox[ "nu" ] )
      }

   #######################   nested CES with 4 inputs   ########################
   } else if( nExog == 4 ) { 
      if( !vrs ) {
         par <- c( par, nu = 1 )
      }

      # derivatives with respect to gamma
      result[ , "gamma" ] <-
         cesDerivCoefN4Gamma( par = par, data = data, xNames = xNames )

      # derivatives with respect to delta_1
      result[ , "delta_1" ] <-
         cesDerivCoefN4Delta1( par = par, data = data, xNames = xNames )
         
      # derivatives with respect to delta_2
      result[ , "delta_2" ] <-
         cesDerivCoefN4Delta2( par = par, data = data, xNames = xNames )

      # derivatives with respect to delta_3
      result[ , "delta_3" ] <- 
         cesDerivCoefN4Delta3( par = par, data = data, xNames = xNames )

      # derivatives with respect to rho_1 and rho_2
      if( returnRho1 ) {
         result[ , "rho_1" ] <-
            cesDerivCoefN4Rho1( par = par, data = data, xNames = xNames )
      }
      if( returnRho2 ) {
         result[ , "rho_2" ] <-
          cesDerivCoefN4Rho2( par = par, data = data, xNames = xNames )
      }

      # derivatives with respect to rho
      if( returnRho ) {
         result[ , "rho" ] <-
            cesDerivCoefN4Rho( par = par, data = data, xNames = xNames )
      }

      # derivatives with respect to nu
      if( vrs ) {
         result[ , "nu" ] <-
            cesDerivCoefN4Nu( par = par, data = data, xNames = xNames )
      }
   } else {
      stop( "the derivatives of the nested CES can be calculated",
         " only for three and four inputs" )
   }

   return( result )
}
