aidsJacobian <- function( allCoef, priceNames, totExpName, data = NULL,
      omitLast = TRUE, alpha0 = 0 ) {
   nObs <- nrow( data )
   coef <- aidsCoef( allCoef )
   nGoods <- length( coef$alpha )
   hom <- all.equal( rowSums( coef$gamma ), rep( 0, nGoods ) ) == TRUE
   sym <- all.equal( coef$gamma, t( coef$gamma ) ) == TRUE
   lnp <- aidsPx( "TL", priceNames, coef = coef, data = data, alpha0 = alpha0 )
   result <- matrix( 0, nrow = nObs * ( nGoods - 1 ),
      ncol = ( nGoods + 2 ) * ( nGoods - 1 ) )
   for( eq in 1:( nGoods - 1 ) ) {
      myRows <- ( ( eq - 1 ) * nObs + 1 ):( eq * nObs )
      # derivatives of alphas
      for( i in 1:( nGoods - 1 ) ) {
         myCol <- ( i - 1 ) * ( nGoods + 2 ) + 1
         result[ myRows, myCol ] <- ( i == eq ) -
            coef$beta[ eq ] *
            ( log( data[[ priceNames[ i ] ]] ) -
            log( data[[ priceNames[ nGoods ] ]] ) )
      }
      # derivatives of betas
      myCol <- ( eq - 1 ) * ( nGoods + 2 ) + 2
      result[ myRows, myCol ] <- log( data[[ totExpName ]] ) - lnp
      # derivatives of gammas
      for( i in 1:( nGoods - 1 ) ) {
         for( j in 1:( nGoods - hom ) ) {
            myCol <- ( i - 1 ) * ( nGoods + 2 ) + 2 + j
            result[ myRows, myCol ] <-
               ( i == eq ) * ( log( data[[ priceNames[ j ] ]] ) -
                  hom * log( data[[ priceNames[ nGoods ] ]] ) ) -
               0.5 * coef$beta[ eq ] *
               ( log( data[[ priceNames[ i ] ]] ) -
                  log( data[[ priceNames[ nGoods ] ]] ) ) *
               ( log( data[[ priceNames[ j ] ]] ) -
                  hom * log( data[[ priceNames[ nGoods ] ]] ) )
         }
      }
   }
   delCols <- NULL
   for( i in 1:( nGoods - 1 ) ) {
      if( hom ) {
         delCols <- c( delCols, i * ( nGoods + 2 ) )
      }
      if( sym && i >= 2 ) {
         for( j in 1:( i - 1 ) ) {
            delCol <- ( i - 1 ) * ( nGoods + 2 ) + 2 + j
            stayCol <- ( j - 1 ) * ( nGoods + 2 ) + 2 + i
            result[ , stayCol ] <- result[ , stayCol ] + result[ , delCol ]
            delCols <- c( delCols, delCol )
         }
      }
   }
   if( !is.null( delCols ) ) {
      result <- result[ , -delCols ]
   }
   return( result )
}
