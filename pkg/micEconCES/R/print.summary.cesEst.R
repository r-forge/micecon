print.summary.cesEst <- function( x, digits = max( 3, getOption( "digits" ) - 3 ),
      ... ) {

   cat( "Estimated CES function with " )
   if( "phi" %in% rownames( coef( x ) ) ){
      cat( "variable " )
   } else {
      cat( "constant " )
   }
   cat( "returns to scale\n\n" )

   cat( "Call:\n" )
   print( x$call )
   cat( "\n" )

   cat( "Estimation by " )
   if( x$method == "Kmenta" ) {
      cat( "the linear Kmenta approximation\n" )
   } else {
      cat( "non-linear least-squares using the '", x$method, "' optimizer\n",
         sep = "" )
      if( ! x$method %in% c( "SANN", "DE" ) ) {
         cat( "Convergence ", ifelse( x$convergence, "", "NOT " ),
            "achieved after ", sep = "" )
         if( length( x$iter ) == 1 ) {
            cat( x$iter, "iterations\n" )
         } else {
            cat( paste( x$iter, names( x$iter ), collapse = " and " ),
               "calls\n" )
         }
      }
      if( !is.null( x$message ) ) {
         cat( "Message:", x$message, "\n" )
      }
   }
   cat( "\n" )

   cat( "Coefficients:\n" )
   printCoefmat( coef( x ), digits = digits )
   cat( "\n" )

   cat( "Residual standard error:", x$sigma, "\n" )
   cat( "Multiple R-squared:", x$r.squared, "\n" )
   cat( "\n" )

   invisible( x )
}
