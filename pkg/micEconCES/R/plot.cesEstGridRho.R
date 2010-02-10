plot.cesEstGridRho <- function( x, ... ) {

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