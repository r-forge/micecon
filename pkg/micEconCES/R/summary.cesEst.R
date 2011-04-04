summary.cesEst <- function( object, rSquaredLog = object$multErr, ... ) {

   # number of observations
   nObs <- length( residuals( object ) )

   # square root of the estimated variance of the random error
   object$sigma <- sqrt( object$rss / nObs )

   # R-squared value
   if( rSquaredLog ) {
      object$r.squared <- rSquared( 
         y = log( fitted( object ) ) + residuals( object ),
         resid = residuals( object ) )
   } else {
      object$r.squared <- rSquared( y = fitted( object ) + residuals( object ),
         resid = residuals( object ) )
   }

   # covariance matrix of the estimated coefficients/parameters
   if( is.null( object$vcov ) ) {
      object$vcov <- object$sigma^2 * object$cov.unscaled
   }

   object$coefficients <- coefTable( coef( object ),
      diag( object$vcov )^0.5, df = Inf )

   class( object ) <- c( "summary.cesEst", class( object ) )
   return( object )
}
