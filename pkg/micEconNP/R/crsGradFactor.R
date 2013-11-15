crsGradFactor <- function( x, varName, onlyOwnLevels = TRUE ) {

   xDat <- x$xz
   yDat <- fitted(x) + residuals(x)

   if( ! varName %in% x$znames ) {
      stop( "the model specified by argument 'x'",
         " does not include an explanatory factor variable '", varName, "'" )
   }
   if( !is.factor( xDat[[ varName ]] ) ) {
      stop( "variable '", varName, "' is not a factor variable" )
   }
   if( ! varName %in% names( xDat ) ) {
      stop( "internal error: the data set in x$eval",
         " does not include a variable '", varName, "'" )
   }

   allLevels <- levels( xDat[[ varName ]] )
   nLevels <- length( allLevels )
   if( nLevels < 2 ) {
      stop( "the factor variable '", varName, "' has must have at least two levels" )
   }
   gradFac <- matrix( NA, nrow = nrow( xDat ), ncol = nLevels - 1 )
   if( is.ordered( xDat[[ varName ]] ) ) {
      colnames( gradFac ) <- paste( varName, ": ", 
         allLevels[ - length( allLevels ) ], " -> ", allLevels[-1], sep = "" )
   } else {
      colnames( gradFac ) <- paste( varName, ": ", allLevels[1], " -> ",
         allLevels[-1], sep = "" )
   }
   for( i in 2:nLevels ) {
      xDatBase <- xDatNew <- xDat
      xDatBase[[ varName ]][] <- levels( xDat[[ varName ]] )[
         ifelse( is.ordered( xDat[[ varName ]] ), i - 1, 1 ) ]
      xDatNew[[ varName ]][] <- levels( xDat[[ varName ]] )[i]
      yHatBase <- predict( x,  newdata = xDatBase )
      yHatNew <- predict( x,  newdata = xDatNew )
      gradFac[ , i - 1 ] <- yHatNew - yHatBase
      if( onlyOwnLevels ) {
         obsOwnLevel <-  xDat[[ varName ]] == allLevels[ i ]
         gradFac[ !obsOwnLevel, i - 1 ] <- NA
      }
   }
   return( gradFac )
}
