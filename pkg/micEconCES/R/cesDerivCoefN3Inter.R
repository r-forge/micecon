cesDerivCoefN3Inter <- function( funcName, par, xNames, data, rhoApprox ) {

      # interpolation if rho and/or rho_1 are close to zero
      rho1 <- par[ "rho_1" ]
      rho <- par[ "rho" ]

      coef00 <- par
      coef01 <- par
      coef10 <- par
      coef11 <- par
      wa <- 1
      wb <- 1
      if( rho == 0 ) {
         wa <- 0
      } else if( abs( rho ) <= rhoApprox ) {
         coef00[ "rho" ] <- coef01[ "rho" ] <- 0
         coef10[ "rho" ] <- coef11[ "rho" ] <- 
            rhoApprox * (-1)^( rho < 0 )
         wa <- abs( rho ) / rhoApprox
      }
      if( rho1 == 0 ) {
         wb <- 0
      } else if( abs( rho1 ) <= rhoApprox ) {
         coef00[ "rho_1" ] <- coef10[ "rho_1" ] <- 0
         coef01[ "rho_1" ] <- coef11[ "rho_1" ] <- 
            rhoApprox * (-1)^( rho1 < 0 )
         wb <- abs( rho1 ) / rhoApprox
      }
      if( wa != 1 && wb != 1 ) {
         result00 <- do.call( funcName, 
            args = list( coef = coef00, data = data, xNames = xNames ) )
      } else {
         result00 <- 0
      }
      if( wa != 1 && wb != 0 ) {
         result01 <- do.call( funcName, 
            args = list( coef = coef01, data = data, xNames = xNames ) )
      } else {
         result01 <- 0
      }
      if( wa != 0 && wb != 1 ) {
         result10 <- do.call( funcName, 
            args = list( coef = coef10, data = data, xNames = xNames ) )
      } else {
         result10 <- 0
      }
      if( wa != 0 && wb != 0 ) {
         result11 <- do.call( funcName, 
            args = list( coef = coef11, data = data, xNames = xNames ) )
      } else {
         result11 <- 0
      }
      result <- 
         wa * wb * result11 + ( 1 - wa ) * wb * result01 +
         wa * ( 1 - wb ) * result10 + ( 1- wa ) * ( 1 - wb ) * result00

   return( result )
}
