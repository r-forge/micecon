plot.cesEst <- function( x, ... ) {

   if( is.null( x$allRhoSum ) ) {
      stop( "the 'plot' method for objects of class 'cesEst' can be applied",
         " only if the CES function was estimated by grid search for 'rho',",
         " i.e. 'cesEst' was called with argument 'rho' set to a vector",
         " with more than one element" )
   }

   rhoVal <- x$allRhoSum$rho
   rssVal <- x$allRhoSum$rss

   plot.default( x = rhoVal, y = rssVal, type = "o", pch = 19,
      xlab = "rho", ylab = "rss", ... )

   # mark estimations that did not converge
   nonConv <- !is.na( x$allRhoSum$convergence )& !x$allRhoSum$convergence
   if( any( nonConv ) ) {
      points( x = rhoVal[ nonConv ], y = rssVal[ nonConv ],
         col = "red", pch = 19 )
   }

   invisible( )
}