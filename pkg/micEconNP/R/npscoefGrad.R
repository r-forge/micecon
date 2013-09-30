npscoefGrad <- function( x, eps = 1e-3, ... ) {
   
   yDat <- fitted(x) + residuals(x)
   xDat <- x$eval$exdat
   zDat <- x$eval$ezdat
   
   nobs <- length( yDat )
   if( nobs != nrow( xDat ) || nobs != nrow( zDat ) ) {
      stop( "internal error: number of observations is not the same for y, x, and z" )
   }
   nx <- ncol( xDat )
   nz <- ncol( zDat )
   
   result <- array( NA, c( nobs, nx + 2, nz ) )
   
   for( i in 1:nz ) {
      zLower <- zUpper <- zDat
      if( is.numeric( zDat[ , i ] ) ) {
         zLower[ , i ] <- zLower[ , i ] - eps / 2
         zUpper[ , i ] <- zUpper[ , i ] + eps / 2
      }
      modelLower <- npscoef( bws = x$bws, 
         tydat = yDat, txdat = xDat, tzdat = zDat, exdat = xDat, ezdat = zLower, betas = TRUE )
      modelUpper <- npscoef( bws = x$bws, 
         tydat = yDat, txdat = xDat, tzdat = zDat, exdat = xDat, ezdat = zUpper, betas = TRUE )
      result[ , 1:(nx+1), i ] <- ( coef( modelUpper ) - coef( modelLower ) ) / eps
      result[ , nx + 2, i ] <- ( fitted( modelUpper ) - fitted( modelLower ) ) / eps
   }
   
   return( result )
}