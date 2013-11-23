npregCv <- function( x ) {
   
   xDat <- x$eval
   yDat <- fitted(x) + residuals(x)
   
   res <- rep( NA, nrow( xDat ) )
   for( i in 1:nrow( xDat ) ) {
      tmp <- npreg( bws = x$bws, tydat = yDat[-i], txdat = xDat[-i,], 
         exdat = xDat[i,] )
      res[ i ] <- yDat[ i ] - fitted( tmp )
   }
   cv <- mean( res^2 )
   attr( cv, "err" ) <- res
   return( cv )
}
