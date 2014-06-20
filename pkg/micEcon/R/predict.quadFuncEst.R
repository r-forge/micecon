predict.quadFuncEst <- function( object, newdata = NULL, ... ) {
   
   if( is.null( newdata ) ) {
      result <- predict( object$est )
   } else {
      result <- quadFuncCalc( xNames = object$xNames, data = newdata,
         coef = coef( object ), shifterNames = object$shifterNames,
         homWeights = object$homWeights )
   }
   return( result )
}
