print.cesEst <- function( x, digits = max( 3, getOption( "digits" ) - 3 ),
      ... ) {

   cat( "Estimated CES function\n\n" )

   cat( "Call:\n" )
   print( x$call )
   cat( "\n" )

   cat( "Coefficients:\n" )
   print( coef( x ), digits = digits )
   cat( "\n" )

   invisible( x )
}