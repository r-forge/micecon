summary.cesEst <- function( object, ... ) {

   # number of observations
   nObs <- length( residuals( object ) )

   # square root of the estimated variance of the random error
   object$sigma <- sqrt( sum( residuals( object )^2 ) / nObs )

   object$coefficients <- coefTable( coef( object ),
      diag( vcov( object ) )^0.5, df = Inf )

   class( object ) <- c( "summary.cesEst", class( object ) )
   return( object )
}