cesDerivCoef <- function( par, data ) {

   # check for VRS
   vrs <- length( par ) >= 4

   # names of coefficients
   coefNames <- cesCoefNames( 2, vrs )

   # derivatives of the CES with respect to the coefficients/parameters
   result <- matrix( NA, nrow = nrow( data ), ncol = length( par ) )
   colnames( result ) <- coefNames
   names( par ) <- coefNames

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
      ( delta * data$x1^(-rho) + ( 1 - delta ) * data$x2^(-rho) )^( -phi / rho )

   # derivatives with respect to delta
   result[ , "delta" ] <- gamma * ( -phi / rho ) *
      ( delta * data$x1^(-rho) + ( 1 - delta ) * data$x2^(-rho) )^( -phi / rho - 1 ) *
      ( data$x1^(-rho) - data$x2^(-rho) )

   # derivatives with respect to rho
   result[ , "rho" ] <- gamma *
      log( delta * data$x1^(-rho) + ( 1 - delta ) * data$x2^(-rho) ) *
      ( delta * data$x1^(-rho) + ( 1 - delta ) * data$x2^(-rho) )^( -phi / rho ) *
      ( phi / rho^2 ) +
      gamma * ( phi / rho ) *
      ( delta * data$x1^(-rho) + ( 1 - delta ) * data$x2^(-rho) )^( -phi / rho - 1 ) *
      ( delta * log( data$x1 ) * data$x1^(-rho) +
         ( 1 - delta ) * log( data$x2 ) * data$x2^(-rho) )

   # derivatives with respect to phi
   if( vrs ) {
      result[ , "phi" ] <- gamma *
         log( delta * data$x1^(-rho) + ( 1 - delta ) * data$x2^(-rho) ) *
         ( delta * data$x1^(-rho) + ( 1 - delta ) * data$x2^(-rho) )^( -phi / rho ) *
         ( -1 / rho )
   }

   return( result )
}
