npplregCv <- function( x ) {
   
   xDat <- x$evalx
   zDat <- x$evalz
   yDat <- fitted(x) + residuals(x)
   
   res <- rep( NA, nrow( xDat ) )
   for( i in 1:nrow( xDat ) ) {
      tmp <- npplreg( bws = x$bw, tydat = yDat[-i], txdat = xDat[-i,], 
         tzdat = zDat[-i,], exdat = xDat[i,], ezdat = zDat[i,] )
      res[ i ] <- yDat[ i ] - fitted( tmp )
   }
   cv <- mean( res^2 )
   attr( cv, "err" ) <- res
   return( cv )
}
