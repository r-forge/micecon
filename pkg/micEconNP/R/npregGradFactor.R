npregGradFactor <- function( x, varName ) {
   data <- x$eval
   grad <- gradients( x )
   if( is.null( grad ) ) {
      stop( "the model specified by argument 'x' does not include gradients" )
   }
   if( ! varName %in% x$xnames ) {
      stop( "the model specified by argument 'x'",
         " does not include an explanatory variable '", varName, "'" )
   }
   if( !is.factor( data[[ varName ]] ) ) {
      stop( "variable '", varName, "' is not a factor variable" )
   }
   if( ! varName %in% names( data ) ) {
      stop( "internal error: the data set in x$eval",
         " does not include a variable '", varName, "'" )
   }
   if( nrow( grad ) != length( data[[ varName ]] ) ) {
      stop( "internal error: gradients have a different number of observations",
         " than variable '", varName, "'" )
   }
   colnames( grad ) <- x$xnames
   allLevels <- levels( data[[ varName ]] )
   nLevels <- length( allLevels )
   if( nLevels < 2 ) {
      stop( "the factor variable '", varName, "' has must have at least two levels" )
   }
   gradFac <- matrix( NA, nrow = nrow( grad ), ncol = nLevels - 1 )
   colnames( gradFac ) <- paste( varName, ": ", allLevels[1], " -> ",
      allLevels[-1], sep = "" )
   for( i in 2:nLevels ) {
      obsSelect <-  data[[ varName ]] == allLevels[ i ]
      gradFac[ obsSelect, i - 1 ] <- grad[ obsSelect, varName ]
   }
   return( gradFac )
}
