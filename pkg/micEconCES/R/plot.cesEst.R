plot.cesEst <- function( x, ... ) {

   if( is.null( x$otherSum ) ) {
      stop( "the 'plot' method for objects of class 'cesEst' can be applied",
         " only if the CES function was estimated by grid search for 'rho',",
         " i.e. 'cesEst' was called with argument 'rho' set to a vector",
         " with more than one element" )
   }

   rhoVal <- x$otherSum$rho
   rssVal <- x$otherSum$rss

   plot.default( x = rhoVal, y = rssVal, type = "o", pch = 19,
      xlab = "rho", ylab = "rss", ... )

   # mark estimations that did not converge
   nonConv <- !is.na( x$otherSum$convergence )& !x$otherSum$convergence
   if( any( nonConv ) ) {
      points( x = rhoVal[ nonConv ], y = rssVal[ nonConv ],
         col = "red", pch = 19 )
   }

   invisible( )
}