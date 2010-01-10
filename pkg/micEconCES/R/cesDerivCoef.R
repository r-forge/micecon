cesDerivCoef <- function( par, data ) {

   # derivatives of the CES with respect to the coefficients/parameters
   result <- matrix( NA, nrow = nrow( data ), ncol = length( par ) )
   colnames( result ) <- c( "gamma", "alpha", "rho", "phi" )[ 1:length( par ) ]

   gamma <- par[ "gamma" ]
   alpha <- par[ "alpha" ]
   rho <- par[ "rho" ]
   if( "phi" %in% names( par ) ) {
      phi <- par[ "phi" ]
   } else {
      phi <- 1
   }

   # derivatives with respect to gamma
   result[ , "gamma" ] <-
      ( alpha * data$x1^(-rho) + ( 1 - alpha ) * data$x2^(-rho) )^( -phi / rho )

   # derivatives with respect to alpha
   result[ , "alpha" ] <- gamma * ( -phi / rho ) *
      ( alpha * data$x1^(-rho) + ( 1 - alpha ) * data$x2^(-rho) )^( -phi / rho - 1 ) *
      ( data$x1^(-rho) - data$x2^(-rho) )

   # derivatives with respect to rho
   result[ , "rho" ] <- gamma *
      log( alpha * data$x1^(-rho) + ( 1 - alpha ) * data$x2^(-rho) ) *
      ( alpha * data$x1^(-rho) + ( 1 - alpha ) * data$x2^(-rho) )^( -phi / rho ) *
      ( phi / rho^2 ) +
      gamma * ( phi / rho ) *
      ( alpha * data$x1^(-rho) + ( 1 - alpha ) * data$x2^(-rho) )^( -phi / rho - 1 ) *
      ( alpha * log( data$x1 ) * data$x1^(-rho) +
         ( 1 - alpha ) * log( data$x2 ) * data$x2^(-rho) )

   # derivatives with respect to phi
   if( "phi" %in% names( par ) ) {
      result[ , "phi" ] <- gamma *
         log( alpha * data$x1^(-rho) + ( 1 - alpha ) * data$x2^(-rho) ) *
         ( alpha * data$x1^(-rho) + ( 1 - alpha ) * data$x2^(-rho) )^( -phi / rho ) *
         ( -1 / rho )
   }

   return( result )
}
