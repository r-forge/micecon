cesRss <- function( par, data ) {
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
   return( sum( ( data$y - yHat )^2 ) )
}
