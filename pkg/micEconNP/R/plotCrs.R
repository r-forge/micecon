plotCrs <- function( x ) {
   
   if( !inherits( x, "crs" ) ) {
      stop( "argument 'x' must be of class 'crs'" )
   }
   
   xVarNames <- colnames( x$deriv.mat )
   
   dat <- eval( x$call$data )
   
   medianDat <- data.frame( no = 1 )

   for( varName in xVarNames ) {
      if( is.factor( dat[[ varName ]] ) ) {
         medianDat[[ varName ]] <-
            names( sort( table( dat[[ varName ]] ), decreasing = TRUE ) )[1]
      } else {
         medianDat[[ varName ]] <- median( dat[[ varName ]] )
      }
   }

   for( varName in xVarNames ) {
      # create data set for making the predictions
      if( is.factor( dat[[ varName ]] ) ) {
         xValAll <- levels( dat[[ varName ]] )
      } else if( is.integer( dat[[ varName ]] ) ) {
         xValAll <- seq( from = min( dat[[ varName ]] ),
            to = max( dat[[ varName ]] ), by = 1 )
      } else {
         xValAll <- seq( from = min( dat[[ varName ]] ),
            to = max( dat[[ varName ]] ), length.out = 100 )
      }
      simDat <- NULL
      for( i in 1:length( xValAll ) ) {
         simDat <- rbind( simDat, medianDat )
      }
      simDat[[ varName ]] <- xValAll
      
      # calculate predicted values
      predVal <- predict( x, newdata = simDat )
      
      # plot
      if( is.factor( dat[[ varName ]] ) ) {
         predVal <- c( predVal )
         names( predVal ) <- simDat[[ varName ]]
         barplot( predVal, main = varName, ylim = range( predVal ) )
      } else {
         plot( simDat[[ varName ]], predVal, main = varName )
      }
   }
   invisible()
}  
