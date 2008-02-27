print.coef.aidsEst <- function( x, ... ) {
   cat( "alpha\n" )
   print( x$alpha )
   cat( "beta\n" )
   print( x$beta )
   cat( "gamma\n" )
   print( x$gamma )
   if( !is.null( x$delta ) ){
      cat( "delta\n" )
      print( x$delta )
   }
   invisible( x )
}
