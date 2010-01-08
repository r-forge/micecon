cesRssDeriv <- function( par, data ) {
   result <- par
   gamma <- par[ "gamma" ]
   alpha <- par[ "alpha" ]
   rho <- par[ "rho" ]
   if( "phi" %in% names( par ) ) {
      phi <- par[ "phi" ]
   } else {
      phi <- 1
   }
   yHat <- gamma *
      ( alpha * data$x1^(-rho) + ( 1 - alpha ) * data$x2^(-rho) )^( -phi / rho )
   resid <- data$y - yHat
   result[ "gamma" ] <- sum( - 2 * resid *
      ( ( alpha * data$x1^(-rho) + ( 1 - alpha ) * data$x2^(-rho) )^( -phi / rho ) ) )
   result[ "alpha" ] <- sum( -2 * resid * ( gamma * ( -phi / rho ) *
      ( alpha * data$x1^(-rho) + ( 1 - alpha ) * data$x2^(-rho) )^( -phi / rho - 1 ) *
      ( data$x1^(-rho) - data$x2^(-rho) ) ) )
   result[ "rho" ] <- sum( -2 * resid * ( gamma *
      log( alpha * data$x1^(-rho) + ( 1 - alpha ) * data$x2^(-rho) ) *
      ( alpha * data$x1^(-rho) + ( 1 - alpha ) * data$x2^(-rho) )^( -phi / rho ) *
      ( phi / rho^2 ) +
      gamma * ( phi / rho ) *
      ( alpha * data$x1^(-rho) + ( 1 - alpha ) * data$x2^(-rho) )^( -phi / rho - 1 ) *
      ( alpha * log( data$x1 ) * data$x1^(-rho) +
         ( 1 - alpha ) * log( data$x2 ) * data$x2^(-rho) ) ) )
   if( "phi" %in% names( par ) ) {
      result[ "phi" ] <- sum( -2 * resid * ( gamma *
         log( alpha * data$x1^(-rho) + ( 1 - alpha ) * data$x2^(-rho) ) *
         ( alpha * data$x1^(-rho) + ( 1 - alpha ) * data$x2^(-rho) )^( -phi / rho ) *
         ( -1 / rho ) ) )
   }

   return( result )
}
