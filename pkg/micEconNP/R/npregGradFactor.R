npregGradFactor <- function( x, varName, onlyOwnLevels = TRUE ) {

   xDat <- x$eval
   yDat <- fitted(x) + residuals(x)
   grad <- gradients( x )
   if( is.null( grad ) ) {
      stop( "the model specified by argument 'x' does not include gradients" )
   }
   if( ! varName %in% x$xnames ) {
      stop( "the model specified by argument 'x'",
         " does not include an explanatory variable '", varName, "'" )
   }
   if( !is.factor( xDat[[ varName ]] ) ) {
      stop( "variable '", varName, "' is not a factor variable" )
   }
   if( ! varName %in% names( xDat ) ) {
      stop( "internal error: the data set in x$eval",
         " does not include a variable '", varName, "'" )
   }
   if( nrow( grad ) != length( xDat[[ varName ]] ) ) {
      stop( "internal error: gradients have a different number of observations",
         " than variable '", varName, "'" )
   }
   colnames( grad ) <- x$xnames
   allLevels <- levels( xDat[[ varName ]] )
   nLevels <- length( allLevels )
   if( nLevels < 2 ) {
      stop( "the factor variable '", varName, "' has must have at least two levels" )
   }
   gradFac <- matrix( NA, nrow = nrow( grad ), ncol = nLevels - 1 )
   colnames( gradFac ) <- paste( varName, ": ", allLevels[1], " -> ",
      allLevels[-1], sep = "" )
   for( i in 2:nLevels ) {
      if( onlyOwnLevels ) {
         obsSelect <-  xDat[[ varName ]] == allLevels[ i ]
         gradFac[ obsSelect, i - 1 ] <- grad[ obsSelect, varName ]
      } else {
         xDatBase <- xDatNew <- xDat
         xDatBase[[ varName ]][] <- levels( xDat[[ varName ]] )[1]
         xDatNew[[ varName ]][] <- levels( xDat[[ varName ]] )[i]
         modelBase <- npreg( bws = x$bws, tydat = yDat, txdat = xDat, 
            exdat = xDatBase )
         modelNew <- npreg( bws = x$bws, tydat = yDat, txdat = xDat,
            exdat = xDatNew )
         gradFac[ , i - 1 ] <- fitted( modelNew ) - fitted( modelBase )         
      }
   }
   return( gradFac )
}
