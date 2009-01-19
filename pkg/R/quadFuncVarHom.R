.quadFuncVarHom <- function( data, xName, homWeights, deflator ) {

   if( is.null( homWeights ) | ! xName %in% names( homWeights ) ) {
      result <- data[[ xName ]]
   } else {
      xOmit <- names( homWeights )[ length( homWeights ) ]
      result <- ( data[[ xName ]] - data[[ xOmit ]] ) / deflator
   }

   return( result )
}