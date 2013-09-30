npregGradFactor <- function( x, varName, data ) {
   grad <- gradients( x )
   if( is.null( grad ) ) {
      stop( "the model specified by argument 'x' does not include gradients" )
   }
   if( ! varName %in% names( data ) ) {
      stop( "the data set specified by argument 'data'",
         " does not include a variable '", varName, "'" )
   }
   if( ! varName %in% x$xnames ) {
      stop( "the model specified by argument 'x'",
         " does not include an explanatory variable '", varName, "'" )
   }
   if( !is.factor( data[[ varName ]] ) ) {
      stop( "variable '", varName, "' is not a factor variable" )
   }
   if( nrow( grad ) != length( data[[ varName ]] ) ) {
      stop( "gradients have a different number of observations",
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
