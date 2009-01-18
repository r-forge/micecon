elas.translogEst <- function( object, data = NULL, dataLogged = NULL, ... ) {

   if( is.null( data ) ) {
      data <- eval( object$call$data )
   }

   if( is.null( dataLogged ) ) {
      if( is.null( object$call$dataLogged ) ) {
         dataLogged <- FALSE
      } else {
         dataLogged <- object$call$dataLogged
      }
   }

   result <- translogEla( xNames = eval( object$call$xNames ), 
      data = data, coef = coef( object ), 
      coefCov = vcov( object ), dataLogged = dataLogged )

   return( result )
}