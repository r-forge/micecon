aidsMono <- function( priceNames, totExpName, data, coef,
      priceIndex = "TL", basePrices = NULL, baseShares = NULL ) {

   if( is.null( coef$alpha0 ) && priceIndex == "TL" ) {
      stop( "argument 'coef' must have element 'alpha0'" )
   }

   result <- list()
   nGoods <- length( priceNames )
   nObs <- nrow( data )

   # fitted shares
   fitted <- aidsCalc( priceNames = priceNames, totExpName = totExpName,
      coef = coef, data = data, priceIndex = priceIndex,
      basePrices = basePrices, baseShares = baseShares )


   # testing for monotonicity
   mono <- array( TRUE, c( nObs ) )
   for( t in 1:nObs ) {
      mono[ t ] <- ( min( fitted$shares[ t, ] ) >= 0 )
   }

   result$mPercent <- 100 * sum( mono ) / nObs
   result$monotony <- mono

   class( result ) <- "aidsMono"
   return( result )
}
