cesRssDeriv <- function( par, data ) {
   result <- par
   gamma <- par[ "gamma" ]
   delta <- par[ "delta" ]
   rho <- par[ "rho" ]
   if( "phi" %in% names( par ) ) {
      phi <- par[ "phi" ]
   } else {
      phi <- 1
   }
   yHat <- cesCalc( xNames = names( data )[ -1 ], data = data, coef = par )
   resid <- data$y - yHat
   result[ "gamma" ] <- sum( - 2 * resid *
      ( ( delta * data$x1^(-rho) + ( 1 - delta ) * data$x2^(-rho) )^( -phi / rho ) ) )
   result[ "delta" ] <- sum( -2 * resid * ( gamma * ( -phi / rho ) *
      ( delta * data$x1^(-rho) + ( 1 - delta ) * data$x2^(-rho) )^( -phi / rho - 1 ) *
      ( data$x1^(-rho) - data$x2^(-rho) ) ) )
   result[ "rho" ] <- sum( -2 * resid * ( gamma *
      log( delta * data$x1^(-rho) + ( 1 - delta ) * data$x2^(-rho) ) *
      ( delta * data$x1^(-rho) + ( 1 - delta ) * data$x2^(-rho) )^( -phi / rho ) *
      ( phi / rho^2 ) +
      gamma * ( phi / rho ) *
      ( delta * data$x1^(-rho) + ( 1 - delta ) * data$x2^(-rho) )^( -phi / rho - 1 ) *
      ( delta * log( data$x1 ) * data$x1^(-rho) +
         ( 1 - delta ) * log( data$x2 ) * data$x2^(-rho) ) ) )
   if( "phi" %in% names( par ) ) {
      result[ "phi" ] <- sum( -2 * resid * ( gamma *
         log( delta * data$x1^(-rho) + ( 1 - delta ) * data$x2^(-rho) ) *
         ( delta * data$x1^(-rho) + ( 1 - delta ) * data$x2^(-rho) )^( -phi / rho ) *
         ( -1 / rho ) ) )
   }

   return( result )
}
