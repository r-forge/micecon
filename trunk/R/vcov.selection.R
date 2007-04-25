vcov.selection <- function( object, part = "full", ... ) {

   if( !( part %in% c( "full", "outcome" ) ) ) {
      stop( "argument 'part' must be either 'full' or 'outcome'" )
   }

   if( object$method == "ml" ){
      result <- vcov.maxlik( object, ... )
   } else if( object$method == "2step" ) {
      result <- object$vcov
   }

   if( part == "outcome" ) {
      if( object$tobitType == 2) {
         i <- c( object$param$index$betaO, object$param$index$Mills )
         result <- result[ i, i ]
      } else if( object$tobitType == 5) {
         i <- c(object$param$index$betaO1, object$param$index$Mills1,
            object$param$index$betaO2, object$param$index$Mills2 )
         result <- result[ i, i ]
      }
   }

   return( result )
}
