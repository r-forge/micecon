snqProfitCalc <- function( pNames, fNames, data, weights, coef, form = 0 ) {

   checkNames( c( pNames, fNames ), names( data ) )

   nNetput <- length( pNames )
   nFix    <- length( fNames )
   nObs    <- nrow( data )

   if( !is.list( coef ) ) {
      stop( "Argument 'coef' must be a list containing the coefficients." )
   }
   if( length( coef$alpha ) != nNetput ) {
      stop( "coef$alpha must have as many elements as argument 'pNames'." )
   }
   if( !is.matrix( coef$beta ) ) {
      stop( "coef$beta must be a matrix." )
   }
   if( nrow( coef$beta ) != ncol( coef$beta ) ) {
      stop( "coef$beta must be a _symmetric_ matrix." )
   }
   if( nrow( coef$beta ) != nNetput ) {
      stop( "coef$beta must have as many rows as argument 'pNames' has elements." )
   }

   normPrice <- numeric( nObs )
   for( i in 1:nNetput ) {
      normPrice <- normPrice + data[[ pNames[ i ] ]] * weights[ i ]
   }

   qNetput <- array( 0, c( nObs, nNetput ) )
   for( i in 1:nNetput ) {
      qNetput[ , i ] <- coef$alpha[ i ]
      for( j in 1:nNetput ) {
         qNetput[ , i ] <- qNetput[ , i ] +
            coef$beta[ i, j ] * data[[ pNames[ j ] ]] / normPrice
         for( k in 1:nNetput ) {
            qNetput[ , i ] <- qNetput[ , i ] -
               0.5 * weights[ i ] * coef$beta[ j, k ] * data[[ pNames[ j ] ]] *
               data[[ pNames[ k ] ]] / normPrice^2
         }
      }
      for( j in 1:nFix ) {
         qNetput[ , i ] <- qNetput[ , i ] +
            coef$delta[ i, j ] * data[[ fNames[ j ] ]]
         for( k in 1:nFix ) {
            if( form == 0 ) {
               qNetput[ , i ] <- qNetput[ , i ] + 0.5 * weights[ i ] *
                  coef$gamma[ j, k ] * data[[ fNames[ j ] ]] * data[[ fNames[ k ] ]]
            } else {
               qNetput[ , i ] <- qNetput[ , i ] + 0.5 * weights[ i ] *
                  coef$gamma[ i, j, k ] * data[[ fNames[ j ] ]] * data[[ fNames[ k ] ]]
            }
         }
      }
   }
   # qNetput <- array(1,c(nObs)) %*% t(coef$alpha) + (P%*%coef$beta) /
   #    (normPrice %*% array(1,c(1,nNetput))) -
   #    0.5 * (diag(P %*% coef$beta %*% t(P)) %*% t(weights)) /
   #    ((normPrice^2) %*% array(1,c(1,nNetput))) +
   #    Z %*% t(coef$delta) +
   #    0.5 * diag(Z %*% coef$gamma %*% t(Z)) %*% t(weights)

   profit <- numeric( nObs )
   for( i in 1:nNetput ) {
      profit <- profit + coef$alpha[ i ] * data[[ pNames[ i ] ]]
      for( j in 1:nNetput ) {
         profit <- profit + 0.5 * coef$beta[ i, j ] * data[[ pNames[ i ] ]] *
            data[[ pNames[ j ] ]] / normPrice
      }
   }
   for( i in 1:nNetput ) {
      for( j in 1:nFix ) {
         profit <- profit + coef$delta[ i, j ] * data[[ pNames[ i ] ]] *
            data[[ fNames[ j ] ]]
      }
   }
   if( form == 0 ) {
      for( i in 1:nFix ) {
         for( j in 1:nFix ) {
            profit <- profit + 0.5 * normPrice * coef$gamma[ i, j ] *
               data[[ fNames[ i ] ]] * data[[ fNames[ j ] ]]
         }
      }
   } else {
       for( i in 1:nNetput ) {
         for( j in 1:nFix ) {
            for( k in 1:nFix ) {
               profit <- profit + 0.5 * normPrice * coef$gamma[ i, j, k ] *
                  data[[ pNames[ i ] ]] * data[[ fNames[ j ] ]] *
                  data[[ fNames[ k ] ]]
            }
         }
      }
   }
   result <- as.data.frame( cbind( qNetput, profit ) )
   names( result ) <- c( paste( "X", 1:nNetput, sep = "" ), "profit" )

   return( result )
}
