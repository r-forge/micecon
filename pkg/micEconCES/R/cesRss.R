cesRss <- function( par, data ) {
   gamma <- par[ "gamma" ]
   delta <- par[ "delta" ]
   rho <- par[ "rho" ]
   if( "phi" %in% names( par ) ) {
      phi <- par[ "phi" ]
   } else {
      phi <- 1
   }
   yHat <- gamma *
      ( delta * data$x1^(-rho) + ( 1 - delta ) * data$x2^(-rho) )^( -phi / rho )
   return( sum( ( data$y - yHat )^2 ) )
}
