coef.selection <- function( object, ... ) {

   coefValues <- object$estimate

   addToCoefNames <- function( prefix, index ) {
      if( !is.null( index ) ){
         names( coefValues )[ index ] <-
            paste( prefix, names( coefValues )[ index ], sep = "" )
      }
      return( coefValues )
   }

   if( "tobit2" %in% class( object ) ) {
      coefValues <- addToCoefNames( "S:", object$param$index$betaS )
      coefValues <- addToCoefNames( "O:", object$param$index$betaO )
   } else if( "tobit5" %in% class( object ) ) {
      coefValues <- addToCoefNames( "S:",  object$param$index$betaS )
      coefValues <- addToCoefNames( "O1:", object$param$index$betaO1 )
      coefValues <- addToCoefNames( "O2:", object$param$index$betaO2 )
   }

   return( coefValues )
}