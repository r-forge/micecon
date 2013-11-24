crsCv <- function( x ) {
   
   yDat <- x$y
   
   mframe <- model.frame( x )
   mform <- x$formula
   
   res <- rep( NA, nrow( mframe ) )
   for( i in 1:nrow( mframe ) ) {
      tmp <- crs( mform, degree = x$degree, segments = x$segments,
         lambda = x$lambda, basis = x$basis, cv = "none", 
         data = mframe[ -i, ] )
      res[ i ] <- yDat[ i ] - predict( tmp, newdata = mframe[ i, ] )
   }
   cv <- mean( res^2 )
   attr( cv, "err" ) <- res
   return( cv )
}
