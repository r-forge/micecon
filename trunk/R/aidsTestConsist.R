aidsTestConsist <- function( pNames, wNames, xtName,
      data = NULL, alpha0, alpha, beta, gamma ) {

   if( length( pNames ) != length( wNames ) ) {
      stop( "arguments 'pNames' and 'wNames' must have the same length." )
   }

   result <- list()
   nGoods <- length( pNames )
   nObs <- length( with( data, get( pNames[ 1 ] ) ) )

   xt <- with( data, get( xtName ) )
   prices <- array( NA, c( nObs, nGoods ) )
   shares <- array( NA, c( nObs, nGoods ) )
   for( i in 1: nGoods ) {
      prices[ , i ] <- with( data, get( pNames[ i ] ) )
      shares[ , i ] <- with( data, get( wNames[ i ] ) )
   }
   fitted <- aidsShares( pNames, xtName, data, alpha0, alpha, beta, gamma )

   # testing for monotonicity
   mono <- array( TRUE, c( nObs ) )
   cMatrices <- list()    # testing for concavity
   conc <- array( TRUE, c( nObs ) )

   lnp <- aidsPx( "TL", pNames, data = data,
      alpha0 = alpha0, alpha=alpha, gamma=gamma )

   for( t in 1:nObs ) {
      mono[ t ] <- ( min( fitted$shares[ t, ] ) >= 0 )
      cMatrices[[ t ]] <- gamma + ( beta %*% t( beta ) ) *
         ( log( xt[ t ] ) - lnp[ t ] ) -
         diag( shares[ t, ] ) + shares[ t, ] %*% t( shares[ t, ] )

  #    for( i in 1:nGoods ) {
  #       conc[ t ] <- ( conc[ t ] & cMatrices[[ t ]][ i, i ] <= 0 )
  #    }
  #    for( i in 2:( nGoods - 1 ) ) {
  #       conc[ t ] <- ( conc[ t ] &
  #      ( det( cMatrices[[ t ]][ 1:i, 1:i ] ) * (-1)^i >= 0 ) )
  #    }
      conc[ t ] <- semidefiniteness( cMatrices[[ t ]][ 1:( nGoods - 1),
         1:( nGoods - 1) ] )$negative
   }
   result$mPercent <- 100 * sum( mono ) / nObs
   result$monotony <- mono
   result$cPercent <- 100 * sum( conc ) / nObs
   result$concavity <- conc
   result$cMatrices <- cMatrices
   return( result )
}
