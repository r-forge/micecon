summary.cesEst <- function( object, ... ) {

   # number of observations
   nObs <- length( residuals( object ) )

   # square root of the estimated variance of the random error
   object$sigma <- sqrt( object$rss / nObs )

   # R-squared value
   object$r.squared <- rSquared( y = fitted( object ) + residuals( object ),
      resid = residuals( object ) )

   # covariance matrix of the estimated coefficients/parameters
   if( is.null( object$vcov ) ) {
      object$vcov <- object$sigma^2 * object$cov.unscaled
   }

   object$coefficients <- coefTable( coef( object ),
      diag( object$vcov )^0.5, df = Inf )

   class( object ) <- c( "summary.cesEst", class( object ) )
   return( object )
}