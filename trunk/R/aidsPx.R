aidsPx <- function( px, pNames, wNames = NULL, data = NULL, base = 1,
   coef = NULL, alpha0 = ifelse( is.null( coef$alpha0 ), 0, coef$alpha0 ) ) {

   nGoods <- length( pNames )
   if( !is.null( wNames ) ) {
      if( nGoods != length( wNames ) && px != "TL" ) {
         stop( "wNames must have as many elements as pNames" )
      }
   }
   nObs <- nrow(  data )
   lnp <- array( 0, c( nObs ))
   if(px=="S") {      # Stone index
      for( i in 1:nGoods ) {
         lnp <- lnp + with( data, get( wNames[ i ] ) ) *
            log( with( data, get( pNames[ i ] ) ) )
      }
   } else if(px=="SL") {     # Stone index with lagged shares
      lnp[ 1 ] <- NA
      for( i in 1:nGoods ) {
         lnp[ 2:nObs ] <- lnp[ 2:nObs ] +
            with( data, get( wNames[ i ] ) )[ 1:(nObs-1) ] *
            log( with( data, get( pNames[ i ] ) ) )[ 2:nObs ]
      }
   } else if(px=="P") {      # log-Paasche index
      for( i in 1:nGoods) {
         lnp <- lnp + with( data, get( wNames[ i ] ) ) *
            log( with( data, get( pNames[ i ] ) ) /
            mean( with( data, get( pNames[ i ] ) )[ base ] ) )
      }
   } else if(px=="L") {      # log-Laspeyres index
      for( i in 1:nGoods) {
         lnp <- lnp + mean( with( data, get( wNames[ i ] ) )[ base ] ) *
            log( with( data, get( pNames[ i ] ) ) )
      }
   } else if(px=="T") {      # Tornqvist index
      for( i in 1:nGoods) {
         lnp <- lnp + c( 0.5 * ( with( data, get( wNames[ i ] ) ) +
            mean( with( data, get( wNames[ i ] ) )[ base ] ) *
            matrix( 1, nrow = nObs ) ) *
            log( with( data, get( pNames[ i ] ) ) /
            mean( with( data, get( pNames[ i ] ) )[ base ] ) ) )
      }
   } else if(px=="TL") {      # Translog index
      lnp <- array( alpha0, c( nObs ) )
      for( i in 1:nGoods ) {
         lnp <- lnp + coef$alpha[ i ] * log( with( data, get( pNames[ i ] ) ) )
         for( j in 1:nGoods ) {
            lnp <- lnp + 0.5 * coef$gamma[ i, j ] *
               log( with( data, get( pNames[ i ] ) ) ) *
               log( with( data, get( pNames[ j ] ) ) )
         }
      }
   } else {
      stop( paste( "The argument 'px' (price index) must be either 'S',",
         "'SL', 'P', 'L', 'T' or 'TL'" ) )
   }
   return( lnp )
}
