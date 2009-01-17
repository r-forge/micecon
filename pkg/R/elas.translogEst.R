elas.translogEst <- function( object, data = NULL, ... ) {

   if( is.null( data ) ) {
      data <- eval( object$call$data )
   }

   result <- translogEla( xNames = eval( object$call$xNames ), 
      data = data, coef = coef( object ), 
      coefCov = vcov( object ) )

   return( result )
}