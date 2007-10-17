coef.selection <- function( object, part="full", ... ) {

   if( !( part %in% c( "full", "outcome" ) ) ) {
      stop( "argument 'part' must be either 'full' or 'outcome'" )
   }
   if("maxLik" %in% class(object))
      coefValues <- coef.maxLik(object)
   else
       coefValues <- object$coefficients
   if( part == "outcome" ) {
      coefValues <- coefValues[ object$param$index$outcome ]
   }

   return( coefValues )
}
