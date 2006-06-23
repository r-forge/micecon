aidsCoef <- function( coef, nGoods, nShifter = 0, cov = NULL, df = 1,
      LA = TRUE, priceNames = NULL, shareNames = NULL ) {
   # nGoods <- -0.5 + ( 2.25 + nrow( array( coef ) ) )^0.5
   if( LA ) {
      M <- matrix( 0, nGoods * ( nGoods + 2 ), ( nGoods - 1 ) * ( nGoods + 2 ) )
      for(i in 1:( nGoods - 1 ) ) {
         M[ i,( i - 1 ) * ( nGoods + 2 ) + 1 ]   <-  1   # alpha[i]
         M[ nGoods, ( i - 1 ) * ( nGoods + 2 ) + 1 ]   <- -1   # alpha[nGoods]
         M[ nGoods + i, ( i - 1 ) * ( nGoods + 2 ) + 2] <-  1   # beta[i]
         M[ nGoods + nGoods, ( i - 1 ) * ( nGoods + 2 ) + 2 ] <- -1
            # beta[ nGoods ]
         for( j in 1:nGoods ) {
            M[ 2 * nGoods + ( i - 1 ) * nGoods + j,
               (i - 1 ) * ( nGoods + 2 ) + 2 + j ] <-  1   # gamma[i,j]
            M[ 2 * nGoods + ( nGoods - 1 ) * nGoods + j,
               ( i - 1 ) * ( nGoods + 2 ) + 2 + j ] <- -1   # gamma[nGoods,j]
         }
      }
      all     <- c( M %*% coef )
      all[nGoods]  <- all[nGoods]+1
      names( all ) <- .aidsCoefNamesAll( nGoods, nShifter )
      alpha   <- all[1:nGoods]
      beta    <- all[(nGoods+1):(2*nGoods)]
      gamma   <- t(array(all[(2*nGoods+1):(nGoods*(nGoods+2))],c(nGoods,nGoods)))
      allcov <- NULL
      stat    <- NULL
      if(!is.null(cov)) {
         allcov   <- M %*% cov %*% t(M)
         rownames( allcov ) <- names( all )
         colnames( allcov ) <- names( all )
         stat     <- coefTable( all, sqrt( diag( allcov ) ), df )
      }
   }
   if( !is.null( shareNames ) ) {
      names( alpha ) <- shareNames
      names( beta ) <- shareNames
      rownames( gamma ) <- shareNames
   }
   if( !is.null( priceNames ) ) {
      colnames( gamma ) <- priceNames
   }
   result <- list()
   result$alpha <- alpha
   result$beta  <- beta
   result$gamma <- gamma
   result$all   <- all
   result$allcov <- allcov
   result$stat  <- stat
   return( result )
}
