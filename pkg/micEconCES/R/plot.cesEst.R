plot.cesEst <- function( x, ... ) {

   if( is.null( x$allRhoSum ) ) {
      stop( "the 'plot' method for objects of class 'cesEst' can be applied",
         " only if the CES function was estimated by grid search for 'rho_1'",
         " or 'rho',",
         " i.e. 'cesEst' was called with argument 'rho1' or 'rho' set to a vector",
         " with more than one element" )
   }

   if( length( dim( x$rssArray ) ) == 3 ) {
      # for three-dimensional grid ssearches
      par( mfcol = c( 3, 1 ), mar = c(0,0,1.5,0) )
      for( i in 1:3 ) {
         if( i == 1 ) {
            xValues <- x$rho1Values
            xLabel <- "rho_1"
            yValues <- x$rho2Values
            yLabel <- "rho_2"
            zValues <- x$rssArray[ , , 
               which( x$rhoValues == coef( x )[ "rho" ] ) ]
         } else if( i == 2 ) {
            xValues <- x$rho1Values
            xLabel <- "rho_1"
            yValues <- x$rhoValues
            yLabel <- "rho"
            zValues <- x$rssArray[ ,
               which( x$rho2Values == coef( x )[ "rho_2" ] ), ]
         } else {
            xValues <- x$rho2Values
            xLabel <- "rho_2"
            yValues <- x$rhoValues
            yLabel <- "rho"
            zValues <- x$rssArray[
               which( x$rho1Values == coef( x )[ "rho_1" ] ), , ]
         }

         # Create a function interpolating colors in the range of specified colors
         jet.colors <- colorRampPalette( c( "green", "red" ) ) 
         # Generate the desired number of colors from this palette
         nbcol <- 100
         color <- jet.colors( nbcol )
         # Compute the z-value at the facet centres
         zfacet <- zValues[ -1, -1 ] + 
            zValues[ -1, -ncol( zValues ) ] + 
            zValues[ -nrow( zValues ), -1 ] + 
            zValues[ -nrow( zValues ), - ncol( zValues ) ]
         # Recode facet z-values into color indices
         facetcol <- cut( log( zfacet ), nbcol )
         # plot
         persp( xValues, yValues, -zValues, 
            phi = 50, theta = -45, expand = 0.75, col = color[ facetcol ],
            xlab = xLabel, ylab = yLabel, zlab = "", ticktype = "detailed",
            main = ifelse( i == 1, "negative sums of squared residuals", "" ),
            ... )
      }
   } else if( is.matrix( x$rssArray ) ) {
      # for two-dimensional grid ssearches
      # Obtain details for graph
      if( is.null( x$rho1Values ) ) {
         xValues <- x$rho2Values
         xLabel <- "rho_2"
         yValues <- x$rhoValues
         yLabel <- "rho"
      } else if( is.null( x$rho2Values ) ) {
         xValues <- x$rho1Values
         xLabel <- "rho_1"
         yValues <- x$rhoValues
         yLabel <- "rho"
      } else if( is.null( x$rhoValues ) ) {
         xValues <- x$rho1Values
         xLabel <- "rho_1"
         yValues <- x$rho2Values
         yLabel <- "rho_2"
      }

      # Create a function interpolating colors in the range of specified colors
      jet.colors <- colorRampPalette( c( "green", "red" ) ) 
      # Generate the desired number of colors from this palette
      nbcol <- 100
      color <- jet.colors( nbcol )
      # Compute the z-value at the facet centres
      zfacet <- x$rssArray[ -1, -1 ] + 
         x$rssArray[ -1, -ncol( x$rssArray ) ] + 
         x$rssArray[ -nrow( x$rssArray ), -1 ] + 
         x$rssArray[ -nrow( x$rssArray ), - ncol( x$rssArray ) ]
      # Recode facet z-values into color indices
      facetcol <- cut( log( zfacet ), nbcol )
      # plot
      persp( xValues, yValues, -x$rssArray, 
         phi = 50, theta = -45, expand = 0.75, col = color[ facetcol ],
         xlab = xLabel, ylab = yLabel, zlab = "", ticktype = "detailed",
         main = "negative sum of squared residuals", ... )
   } else if( is.null( x$rssArray ) ) { 
      # for one-dimensional grid searches
      if( !is.null( x$allRhoSum[[ "rho1" ]] ) ) {
         rhoVal <- x$allRhoSum[[ "rho1" ]]
         rhoName <- "rho_1"
      } else if( !is.null( x$allRhoSum[[ "rho2" ]] ) ) {
         rhoVal <- x$allRhoSum[[ "rho2" ]]
         rhoName <- "rho_2"
      } else if( !is.null( x$allRhoSum[[ "rho" ]] ) ) {
         rhoVal <- x$allRhoSum[[ "rho" ]]
         rhoName <- "rho"
      } else {
         stop( "either 'x$allRhoSum$rho1', 'x$allRhoSum$rho2', or",
            " 'x$allRhoSum$rho' must be non-NULL" )
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
