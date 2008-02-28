aidsConcav <- function( priceNames, totExpName, coef, data,
      shareNames = NULL ) {

   if( !is.null( shareNames ) && length( priceNames ) != length( shareNames ) ) {
      stop( "arguments 'priceNames' and 'shareNames' must have the same length" )
   }
   if( is.null( coef$alpha0 ) ) {
      stop( "argument 'coef' must have element 'alpha0'" )
   }

   result <- list()
   nGoods <- length( priceNames )
   nObs <- nrow( data )

   xt <- data[[ totExpName ]]
   shareMat <- array( NA, c( nObs, nGoods ) )
   for( i in 1: nGoods ) {
      if( !is.null( shareNames ) ) {
         shareMat[ , i ] <- data[[ shareNames[ i ] ]]
      }
   }
   fitted <- aidsCalc( priceNames, totExpName, data = data,
      coef = coef )
   if( is.null( shareNames ) ) {
      shareMat <- as.matrix( fitted$shares )
   }

   # checking concavity
   cMatrices <- list()    
   conc <- array( TRUE, c( nObs ) )

   lnp <- aidsPx( "TL", priceNames, data = data, coef = coef )

   for( t in 1:nObs ) {
      cMatrices[[ t ]] <- coef$gamma + ( coef$beta %*% t( coef$beta ) ) *
         ( log( xt[ t ] ) - lnp[ t ] ) -
         diag( shareMat[ t, ] ) + shareMat[ t, ] %*% t( shareMat[ t, ] )

      conc[ t ] <- semidefiniteness( cMatrices[[ t ]][ 1:( nGoods - 1),
         1:( nGoods - 1) ] )$negative
   }
   result$cPercent <- 100 * sum( conc ) / nObs
   result$concavity <- conc
   result$cMatrices <- cMatrices
   class( result ) <- "aidsConcav"
   return( result )
}
