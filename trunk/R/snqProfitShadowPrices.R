## ===== calculation of elasticities from beta matrix ===
snqProfitShadowPrices <- function( coef, price, fix, normPrices ) {
   if( is.null( dim( price ) ) ) {
      price <- array( price, c( 1, length( price ) ) )
   }
   if( is.null( dim( fix ) ) ) {
      fix <- array( fix, c( 1, length( fix ) ) )
   }
   nNetput <- nrow( coef$delta )
   nFix    <- ncol( coef$delta )
   nObs    <- nrow( fix )
   shadowPrices <- array( 0, c( nObs, nFix ) )

   for( i in 1:nFix ) {
      for( j in 1:nNetput ) {
         shadowPrices[ , i ] <- shadowPrices[ , i ] + coef$delta[ j, i ] * price[ , j ]
      }
      for( j in 1:nFix ) {
         shadowPrices[ , i ] <- shadowPrices[ , i ] + normPrices *
            coef$gamma[ j, i ] * fix[ , j ]
      }
   }
   return( shadowPrices )
}
