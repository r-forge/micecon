cesEst <- function( yName, xNames, data, vrs = FALSE, ... ) {

   # y = gamma * ( alpha * x1^rho + ( 1 - alpha ) * x2^rho )^(phi/rho)
   # s = 1 / ( 1 - rho )

   checkNames( c( yName, xNames ), names( data ) )

   if( length( xNames ) != 2 ) {
      stop( "currently, argument 'xNames' must contain exactly",
         " two variable names" )
   }

   # prepare data for estimation
   estData <- data.frame( y = data[[ yName ]],
      x1 = data[[ xNames[ 1 ] ]], x2 = data[[ xNames[ 2 ] ]] )

   # start values
   startVal <- c( gamma = sqrt( mean( estData$y ) ),
      alpha = 0.5, rho = 0.5 )
   if( vrs ) {
      startVal <- c( startVal, phi = 1 )
   }

   cesRSS <- function( par, data ) {
      gamma <- par[ "gamma" ]
      alpha <- par[ "alpha" ]
      rho <- par[ "rho" ]
      if( "phi" %in% names( par ) ) {
         phi <- par[ "phi" ]
      } else {
         phi <- 1
      }
      yHat <- gamma *
         ( alpha * data$x1^rho + ( 1 - alpha ) * data$x2^rho )^( phi / rho )
      return( sum( ( data$y - yHat )^2 ) )
   }

   result <- optim( startVal, cesRSS, data = estData, ... )

   # nonlinear least squares
#    result$nls <- nls(
#       y ~ gamma * ( alpha * x1^rho + ( 1 - alpha ) * x2^rho )^(1/rho),
#       data = estData, start = result$startVal, trace = TRUE,
#       algorithm = "port", lower = c( -Inf, 0.01 , -Inf ),
#       upper = c( Inf, 0.99, Inf ) )

   class( result ) <- c( "cesEst", class( result ) )
   return( result )
}

