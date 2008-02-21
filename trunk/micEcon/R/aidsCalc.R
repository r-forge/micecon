aidsCalc <- function( priceNames, totExpName, data = NULL, priceIndex = "TL", lnp = NULL,
   coef = NULL ) {

   # check argument 'coef' (coefficients)
   if( !is.null( coef ) ){
      coefCheckResult <- .aidsCheckCoef( coef, variables = list(
         list( length( priceNames ), "prices", "goods"  ) ) )
      if( !is.null( coefCheckResult ) ){
         stop( coefCheckResult )
      }
   }

   # check whether the price index is provided if it should not be in translog form
   if( ! priceIndex %in% c( "TL", "S" ) && is.null( lnp ) ) {
      stop( "at the moment only the translog (TL) and Stone (S) price index work",
         " if argument 'lnp' is not specified" )
   }

   # calculate price index if it isn't provided
   if( is.null( lnp ) && priceIndex == "TL" ) {
      if( is.null( coef$alpha0 ) ) {
         stop( "calculations with the translog (TL) price index require",
            " coefficient alpha_0 (coef$alpha0)" )
      }
      lnp <- aidsPx( priceIndex, priceNames, data = data,
         alpha0 = coef$alpha0, coef = coef )
   }

   # number of goods
   nGoods <- length( priceNames )

   shareData <- as.data.frame( matrix( 0, nrow = nrow( data ), ncol = nGoods ) )
   names( shareData ) <- paste( "w", as.character( 1:nGoods ), sep = "" )
   rownames( shareData ) <- rownames( data )
   quant <- as.data.frame( matrix( 0, nrow = nrow( data ), ncol = nGoods ) )
   names( quant ) <- paste( "q", as.character( 1:nGoods ), sep = "" )
   rownames( quant ) <- rownames( data )
   if( !is.null( lnp ) ) {
      for( i in 1:nGoods ) {
         shareData[ , i ] <- coef$alpha[ i ] + coef$beta[ i ] *
            ( log( data[[ totExpName ]] ) - lnp )
         for( j in 1:nGoods ) {
            shareData[ , i ] <- shareData[ , i ] + coef$gamma[ i, j ] *
               log( data[[ priceNames[ j ] ]] )
         }
      }
   } else if( priceIndex == "S" ) {
      for( i in 1:nrow( data ) ) {
         logPrices <- log( as.numeric( data[ i, priceNames ] ) )
         logTotExp <- log( data[ i, totExpName ] )
         shareData[ i, ] <-
            solve( diag( nGoods ) + coef$beta %*% t( logPrices ),
               coef$alpha + coef$gamma %*% logPrices + coef$beta * logTotExp )
      }
   } else {
      stop( "internal error" )
   }
   for( i in 1:nGoods ) {
      quant[ , i ] <- shareData[ , i ] * data[[ totExpName ]] / data[[ priceNames[ i ] ]]
   }
   result <- list()
   result$shares <- shareData
   result$quant  <- quant
   return( result )
}
