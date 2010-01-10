print.summary.cesEst <- function( x, digits = max( 3, getOption( "digits" ) - 3 ),
      ... ) {

   cat( "Estimated CES function with " )
   if( is.null( x$call$vrs ) || !x$call$vrs ){
      cat( "constant " )
   } else {
      cat( "variable " )
   }
   cat( "returns to scale\n\n" )

   cat( "Call:\n" )
   print( x$call )
   cat( "\n" )

   cat( "Estimation method: ", x$method, "\n\n" )

   cat( "Coefficients:\n" )
   printCoefmat( coef( x ), digits = digits )
   cat( "\n" )

   cat( "Residual standard error:", x$sigma, "\n" )
   cat( "\n" )

   invisible( x )
}
