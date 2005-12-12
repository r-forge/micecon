.snqProfitModelData <- function( data, weights, pNames, qNames, fNames,
   ivNames, form, fixedScale ){

   nNetput <- length( qNames )  # number of netputs
   nFix    <- length( fNames )  # number of fixed inputs
   nIV     <- length( ivNames )  # number of fixed inputs
   nObs    <- nrow( data )      # number of observations

   ## price index for normalization
   result <- data.frame( nr = 1:nObs, normPrice = 0 )
   for( i in 1:nNetput ) {
      result$normPrice <- result$normPrice +
         data[[ pNames[ i ] ]] * weights[ i ]
   }

   ## real/normalized netput prices and netput quantities
   for( i in 1:nNetput ) {
      result[[ paste( "pr", as.character( i ), sep = "" ) ]] <-
         data[[ pNames[ i ] ]] / result$normPrice
      result[[ paste( "q", as.character( i ), sep = "" ) ]] <-
         data[[ qNames[ i ] ]]
   }

   ## quadratic netput prices
   for( i in 1:nNetput ) {
      for( j in 1:nNetput ) {
         for( k in 1:nNetput ) {
            result[[ paste( "pq", as.character( i ), ".", as.character( j ), ".",
               as.character( k ), sep = "" ) ]] <-
               -0.5 * weights[ i ] * data[[ pNames[ j ] ]] *
               data[[ pNames[ k ] ]] / result$normPrice^2
         }
      }
   }

   ## quasi-fix inputs
   if( nFix > 0 ) {
      for( i in 1:nFix ) {
         result[[ paste( "f", as.character( i ), sep = "" ) ]] <-
            data[[ fNames[ i ] ]] / fixedScale[ i ]
      }
      ## quadratic quasi-fix inputs
      for( i in 1:nNetput ) {
         for( j in 1:nFix ) {
            for( k in 1:nFix ) {
               result[[ paste( "fq", as.character( i ), ".", as.character( j ), ".",
               as.character( k ), sep = "" ) ]] <-
                  0.5 * ifelse( form == 0, weights[ i ], 1 ) *
                  ( data[[ fNames[ j ] ]] / fixedScale[ j ] ) *
                  ( data[[ fNames[ k ] ]] / fixedScale[ k ] )
            }
         }
      }
   }

   ## instrumental variables
   if( nIV > 0 ) {
      for( i in 1:nIV ) {
         result[[ paste( "iv", as.character( i ), sep = "" ) ]] <-
            data[[ ivNames[ i ] ]] / mean( data[[ ivNames[ i ] ]] )
      }
   }
   return( result )
}
