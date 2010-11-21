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
         coefArray[ "rho", 1, 1 ] <- coefArray[ "rho", 1, 2 ] <- 0
         coefArray[ "rho", 2, 1 ] <- coefArray[ "rho", 2, 2 ] <- 
            rhoApprox * (-1)^( rho < 0 )
         wa <- abs( rho ) / rhoApprox
      }
      if( rho1 == 0 ) {
         wb <- 0
      } else if( abs( rho1 ) <= rhoApprox ) {
         coefArray[ "rho_1", 1, 1 ] <- coefArray[ "rho_1", 2, 1 ] <- 0
         coefArray[ "rho_1", 1, 2 ] <- coefArray[ "rho_1", 2, 2 ] <- 
            rhoApprox * (-1)^( rho1 < 0 )
         wb <- abs( rho1 ) / rhoApprox
      }
      if( wa != 1 && wb != 1 ) {
         result00 <- do.call( funcName, 
            args = list( coef = coefArray[ , 1, 1 ], data = data, xNames = xNames ) )
      } else {
         result00 <- 0
      }
      if( wa != 1 && wb != 0 ) {
         result01 <- do.call( funcName, 
            args = list( coef = coefArray[ , 1, 2 ], data = data, xNames = xNames ) )
      } else {
         result01 <- 0
      }
      if( wa != 0 && wb != 1 ) {
         result10 <- do.call( funcName, 
            args = list( coef = coefArray[ , 2, 1 ], data = data, xNames = xNames ) )
      } else {
         result10 <- 0
      }
      if( wa != 0 && wb != 0 ) {
         result11 <- do.call( funcName, 
            args = list( coef = coefArray[ , 2, 2 ], data = data, xNames = xNames ) )
      } else {
         result11 <- 0
      }
      result <- 
         wa * wb * result11 + ( 1 - wa ) * wb * result01 +
         wa * ( 1 - wb ) * result10 + ( 1- wa ) * ( 1 - wb ) * result00

   return( result )
}
