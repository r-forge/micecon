aidsPx <- function( priceIndex, priceNames, shareNames = NULL, data = NULL, base = 1,
   coef = NULL, shifterNames = NULL ) {

   if( priceIndex == "TL" ){
      if( is.null( coef ) ) {
         stop( "argument 'coef' must be specified to calculate the translog",
            " price index" )
      } else {
         coefCheckResult <- .aidsCheckCoef( coef, variables = list(
            list( length( priceNames ), "priceNames", "goods" ),
            list( ifelse( is.null( shareNames ), NA, length( shareNames ) ), 
               "shareNames", "goods" ) ) )
         if( !is.null( coefCheckResult ) ){
            stop( coefCheckResult )
         }
         if( is.null( coef$alpha0 ) ) {
            stop( "argument 'coef' must have element 'alpha0'" )
         }
      }
   } else {
      if( is.null( shareNames ) ) {
         stop( "argument 'shareNames' must must be specified to calculate",
            " price index '", priceIndex, "'" )
      }
   }

   nGoods <- length( priceNames )
   nShifter <- length( shifterNames )
   nObs <- nrow(  data )
   lnp <- array( 0, c( nObs ))
   if(priceIndex=="S") {      # Stone index
      for( i in 1:nGoods ) {
         lnp <- lnp + data[[ shareNames[ i ] ]] * log( data[[ priceNames[ i ] ]] )
      }
   } else if(priceIndex=="SL") {     # Stone index with lagged shares
      lnp[ 1 ] <- NA
      for( i in 1:nGoods ) {
         lnp[ 2:nObs ] <- lnp[ 2:nObs ] +
            data[[ shareNames[ i ] ]][ 1:(nObs-1) ] *
            log( data[[ priceNames[ i ] ]][ 2:nObs ] )
      }
   } else if(priceIndex=="P") {      # log-Paasche index
      for( i in 1:nGoods) {
         lnp <- lnp + data[[ shareNames[ i ] ]] * log( data[[ priceNames[ i ] ]] /
            mean( data[[ priceNames[ i ] ]][ base ] ) )
      }
   } else if(priceIndex=="L") {      # log-Laspeyres index
      for( i in 1:nGoods) {
         lnp <- lnp + mean( data[[ shareNames[ i ] ]][ base ] ) *
            log( data[[ priceNames[ i ] ]] )
      }
   } else if(priceIndex=="T") {      # Tornqvist index
      for( i in 1:nGoods) {
         lnp <- lnp + c( 0.5 * ( data[[ shareNames[ i ] ]] +
            mean( data[[ shareNames[ i ] ]][ base ] ) *
            matrix( 1, nrow = nObs ) ) * log( data[[ priceNames[ i ] ]] /
            mean( data[[ priceNames[ i ] ]][ base ] ) ) )
      }
   } else if(priceIndex=="TL") {      # Translog index
      lnp <- array( coef$alpha0, c( nObs ) )
      for( i in 1:nGoods ) {
         lnp <- lnp + coef$alpha[ i ] * log( data[[ priceNames[ i ] ]] )
         for( j in 1:nGoods ) {
            lnp <- lnp + 0.5 * coef$gamma[ i, j ] *
               log( data[[ priceNames[ i ] ]] ) *
               log( data[[ priceNames[ j ] ]] )
         }
         if( nShifter > 0 ){
            for( j in 1:nShifter ) {
               lnp <- lnp + coef$delta[ i, j ] * data[[ shifterNames[ j ] ]] *
                  log( data[[ priceNames[ i ] ]] )
            }
         }
      }
   } else {
      stop( "the argument 'priceIndex' (price index) must be either 'S',",
         " 'SL', 'P', 'L', 'T' or 'TL'" )
   }
   return( lnp )
}
