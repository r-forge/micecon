snqProfitShadowPrices <- function( pNames, fNames, data, weights, coef, form = 0 ) {

   checkNames( c( pNames, fNames ), names( data ) )

   nNetput <- length( pNames )
   nFix    <- length( fNames )
   nObs    <- nrow( data )

   snqProfitTestCoef( nNetput, nFix, coef, form = form,
      coefNames = c( "delta", "gamma" ) )

   normPrice <- numeric( nObs )
   for( i in 1:nNetput ) {
      normPrice <- normPrice + data[[ pNames[ i ] ]] * weights[ i ]
   }

   shadowPrices <- array( 0, c( nObs, nFix ) )
   for( i in 1:nFix ) {
      for( j in 1:nNetput ) {
         shadowPrices[ , i ] <- shadowPrices[ , i ] + coef$delta[ j, i ] *
            data[[ pNames[ j ] ]]
      }
      for( j in 1:nFix ) {
         shadowPrices[ , i ] <- shadowPrices[ , i ] + normPrice *
            coef$gamma[ j, i ] * data[[ fNames[ j ] ]]
      }
   }
   result <- as.data.frame( shadowPrices )
   names( result ) <- fNames

   return( result )
}
