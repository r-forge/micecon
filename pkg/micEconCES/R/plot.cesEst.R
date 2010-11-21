plot.cesEst <- function( x, ... ) {

   if( is.null( x$allRhoSum ) ) {
      stop( "the 'plot' method for objects of class 'cesEst' can be applied",
         " only if the CES function was estimated by grid search for 'rho_1'",
         " or 'rho',",
         " i.e. 'cesEst' was called with argument 'rho1' or 'rho' set to a vector",
         " with more than one element" )
   }

   if( !is.null( x$allRhoSum[[ "rho1" ]] ) && 
         !is.null( x$allRhoSum[[ "rho" ]] ) ) {
      # for two-dimensional grid ssearches
      # Create a function interpolating colors in the range of specified colors
      jet.colors <- colorRampPalette( c( "green", "red" ) ) 
      # Generate the desired number of colors from this palette
      nbcol <- 100
      color <- jet.colors( nbcol )
      # Compute the z-value at the facet centres
      zfacet <- x$rssMatrix[ -1, -1 ] + 
         x$rssMatrix[ -1, -ncol( x$rssMatrix ) ] + 
         x$rssMatrix[ -nrow( x$rssMatrix ), -1 ] + 
         x$rssMatrix[ -nrow( x$rssMatrix ), - ncol( x$rssMatrix ) ]
      # Recode facet z-values into color indices
      facetcol <- cut( log( zfacet ), nbcol )
      # plot
      persp( x$rho1Values, x$rhoValues, -x$rssMatrix, 
         phi = 50, theta = -45, expand = 0.75, col = color[ facetcol ],
         xlab = "rho_1", ylab = "rho", zlab = "- RSS", ticktype = "detailed" )
   } else { 
      # for one-dimensional grid searches
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
   }

   invisible( )
}
