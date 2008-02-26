aidsCalc <- function( priceNames, totExpName, coef, data,
      priceIndex = "TL", basePrices = NULL, baseShares = NULL ) {

   # check argument 'coef' (coefficients)
   coefCheckResult <- .aidsCheckCoef( coef, variables = list(
      list( length( priceNames ), "prices", "goods"  ) ) )
   if( !is.null( coefCheckResult ) ){
      stop( coefCheckResult )
   }

   # checking argument 'data'
   if( class( data ) != "data.frame" ) {
      stop( "argument 'data' must be a data frame" )
   }

   # checking (mainly) argument 'priceIndex'
   if( is.character( priceIndex ) ) {
      if( ! priceIndex %in% c( "TL", "S", "Ls" ) ) {
         stop( "at the moment, argument 'priceIndex' must be either",
            " 'TL' (translog), 'S' (Stone), 'Ls' (Laspeyres, simplified), or a numeric vector",
            " providing the log values of the price index" )
      }
      if( priceIndex == "TL" && is.null( coef$alpha0 ) ) {
         stop( "calculations with the translog (TL) price index require",
            " coefficient alpha_0 (coef$alpha0)" )
      }
   } else if( is.numeric( priceIndex ) ) {
      if( length( priceIndex ) != nrow( data ) ) {
         stop( "if argument 'priceIndex' provides the values",
            " of the log price index,",
            " it must have the same length as there are observations",
            " in argument 'data'" )
      }
   } else {
      stop( "argument 'priceIndex' must be either a charater string",
         " or a numeric vector" )
   }

   # tests for arguments basePrices and baseShares
   if( is.character( priceIndex ) ) {
      if( priceIndex == "Ls" ) {
         # basePrices
         if( is.null( basePrices ) ) {
#             stop( "calculations with simplified Laspeyres ('Ls') price index require",
#                " argument 'basePrices'" )
         }
         if( ! is.numeric( basePrices ) ) {
            stop( "argument 'basePrices' must be numeric" )
         }
         if( length( basePrices ) != length( priceNames ) ) {
            stop( "arguments 'basePrices' and 'priceNames' must have",
               " the same length" )
         }
         # baseShares
         if( is.null( baseShares ) ) {
            stop( "calculations with simplified Laspeyres ('Ls') price index require",
               " argument 'baseShares'" )
         }
         if( ! is.numeric( baseShares ) ) {
            stop( "argument 'baseShares' must be numeric" )
         }
         if( length( baseShares ) != length( priceNames ) ) {
            stop( "arguments 'baseShares' and 'priceNames' must have",
               " the same length" )
         }
      }
   }

   if( is.character( priceIndex ) ) {
      if( priceIndex == "TL" ) {
         # calculation of translog price index
         priceIndex <- aidsPx( priceIndex, priceNames, data = data, coef = coef )
      } else if( priceIndex == "Ls" ) {
         # calculation of simplified Laspeyres price index
         priceIndex <- aidsPx( priceIndex, priceNames, data = data,
            coef = coef, base = list( prices = basePrices, shares = baseShares ) )
      }
   }

   # number of goods
   nGoods <- length( priceNames )

   shareData <- as.data.frame( matrix( 0, nrow = nrow( data ), ncol = nGoods ) )
   names( shareData ) <- paste( "w", as.character( 1:nGoods ), sep = "" )
   rownames( shareData ) <- rownames( data )
   quant <- as.data.frame( matrix( 0, nrow = nrow( data ), ncol = nGoods ) )
   names( quant ) <- paste( "q", as.character( 1:nGoods ), sep = "" )
   rownames( quant ) <- rownames( data )
   if( is.numeric( priceIndex ) ) {
      for( i in 1:nGoods ) {
         shareData[ , i ] <- coef$alpha[ i ] + coef$beta[ i ] *
            ( log( data[[ totExpName ]] ) - priceIndex )
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
