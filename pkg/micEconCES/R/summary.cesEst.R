summary.cesEst <- function( object, rSquaredLog = object$multErr, ... ) {

   # number of observations
   nObs <- length( residuals( object ) )

   # square root of the estimated variance of the random error
   object$sigma <- sqrt( object$rss / nObs )

   # R-squared value
   fitVal <- fitted( object ) # the returned fitted values are never logged
   if( rSquaredLog ) {
      if( object$multErr ) {
         yVal <- log( fitVal ) + residuals( object )  # log( y )
      } else {
         yVal <- log( fitVal + residuals( object ) )  # log( y )
      }
      res <- yVal - log( fitVal )
   } else {
      if( object$multErr ) {
         yVal <- fitVal * exp( residuals( object ) )  # y
      } else {
         yVal <- fitVal + residuals( object )         # y
      }
      res <- yVal - fitVal
   }
   object$r.squared <- rSquared( y = yVal, resid = res )

   # covariance matrix of the estimated coefficients/parameters
   if( is.null( object$vcov ) ) {
      object$vcov <- object$sigma^2 * object$cov.unscaled
   }

   object$coefficients <- coefTable( coef( object ),
      diag( object$vcov )^0.5, df = Inf )

   class( object ) <- c( "summary.cesEst", class( object ) )
   return( object )
}
