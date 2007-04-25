coef.summary.selection <- function( object, part="full", ... ) {

   if( !( part %in% c( "full", "outcome" ) ) ) {
      stop( "argument 'part' must be either 'full' or 'outcome'" )
   }

   result <- object$estimate

   if( part == "outcome" ) {
      if( object$tobitType == 2) {
         result <- result[ c( object$param$index$betaO,
            object$param$index$Mills ), ]
      } else if( object$tobitType == 5) {
         result <- result[ c(object$param$index$betaO1,
            object$param$index$Mills1, object$param$index$betaO2,
            object$param$index$Mills2 ), ]
      }
   }

   return( result )
}
