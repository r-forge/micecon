cesInterN3 <- function( funcName, par, xNames, data, rhoApprox ) {

      # interpolation if rho and/or rho_1 are close to zero
      rho1 <- par[ "rho_1" ]
      rho <- par[ "rho" ]

      coefArray <- array( par, c( length( par ), 2, 2 ) )
      dimnames( coefArray ) <- list( names( par ), 
         c( "rho = 0", "rho = E" ),  c( "rho_1 = 0", "rho_1 = E" ) )
      wa <- 1
      wb <- 1
      if( rho == 0 ) {
         wa <- 0
      } else if( abs( rho ) <= rhoApprox ) {
         coefArray[ "rho", "rho = 0", ] <- 0
         coefArray[ "rho", "rho = E", ] <- rhoApprox * (-1)^( rho < 0 )
         wa <- abs( rho ) / rhoApprox
      }
      if( rho1 == 0 ) {
         wb <- 0
      } else if( abs( rho1 ) <= rhoApprox ) {
         coefArray[ "rho_1", , "rho_1 = 0" ] <- 0
         coefArray[ "rho_1", , "rho_1 = E" ] <- rhoApprox * (-1)^( rho1 < 0 )
         wb <- abs( rho1 ) / rhoApprox
      }
      result <- array( 0, c( nrow( data ), 2, 2 ) )
      if( wa != 1 && wb != 1 ) {
         result[ , 1, 1 ] <- do.call( funcName, 
            args = list( coef = coefArray[ , 1, 1 ], data = data, xNames = xNames ) )
      }
      if( wa != 1 && wb != 0 ) {
         result[ , 1, 2 ] <- do.call( funcName, 
            args = list( coef = coefArray[ , 1, 2 ], data = data, xNames = xNames ) )
      }
      if( wa != 0 && wb != 1 ) {
         result[ , 2, 1 ] <- do.call( funcName, 
            args = list( coef = coefArray[ , 2, 1 ], data = data, xNames = xNames ) )
      }
      if( wa != 0 && wb != 0 ) {
         result[ , 2, 2 ] <- do.call( funcName, 
            args = list( coef = coefArray[ , 2, 2 ], data = data, xNames = xNames ) )
      }
      result <- 
         wa * wb * result[ , 2, 2 ] + ( 1 - wa ) * wb * result[ , 1, 2 ] +
         wa * ( 1 - wb ) * result[ , 2, 1 ] + ( 1- wa ) * ( 1 - wb ) * result[ , 1, 1 ]

   return( result )
}
