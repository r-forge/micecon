summary.cesEst <- function( object, ... ) {

   object$coefficients <- coefTable( coef( object ),
      diag( vcov( object ) )^0.5, df = Inf )

   class( object ) <- c( "summary.cesEst", class( object ) )
   return( object )
}