plot.cesEst <- function( x, ... ) {

   if( is.null( x$allRhoSum ) ) {
      stop( "the 'plot' method for objects of class 'cesEst' can be applied",
         " only if the CES function was estimated by grid search for 'rho_1'",
         " or 'rho',",
         " i.e. 'cesEst' was called with argument 'rho1' or 'rho' set to a vector",
         " with more than one element" )
   }

   if( !is.null( x$allRhoSum[[ "rho1" ]] ) ) {
      rhoVal <- x$allRhoSum[[ "rho1" ]]
      rhoName <- "rho_1"
   } else if( !is.null( x$allRhoSum[[ "rho" ]] ) ) {
      rhoVal <- x$allRhoSum[[ "rho" ]]
      rhoName <- "rho"
   } else {
      stop( "either 'x$allRhoSum$rho1' or 'x$allRhoSum$rho' must be non-NULL" )
   }
   rssVal <- x$allRhoSum$rss

   plot.default( x = rhoVal, y = rssVal, type = "o", pch = 19,
      xlab = rhoName, ylab = "rss", ... )

   # mark estimations that did not converge
   nonConv <- !is.na( x$allRhoSum$convergence )& !x$allRhoSum$convergence
   if( any( nonConv ) ) {
      points( x = rhoVal[ nonConv ], y = rssVal[ nonConv ],
         col = "red", pch = 19 )
   }

   invisible( )
}
