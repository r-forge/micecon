aidsElaJacobian <- function( coef, share, price = NULL, formula = "AIDS" ) {

   if( length( coef$alpha ) != length( coef$beta ) ) {
      stop( "arguments 'alpha' and 'beta' must have the same length" )
   } else if( nrow( coef$gamma ) != ncol( coef$gamma ) ) {
      stop( "argument 'gamma' must be a square matrix" )
   } else if( length( coef$alpha ) != nrow( coef$gamma ) ) {
      stop( "number of rows of argument 'gamma' must be equal",
         " to the length of argument 'alpha'" )
   } else if(  length( coef$alpha ) != length( share ) ) {
      stop( "arguments 'alpha' and 'share' must have the same length" )
   } else if(  length( coef$alpha ) != length( price ) && !is.null( price ) ) {
      stop( "arguments 'alpha' and 'price' must have the same length" )
   }

   if( formula %in% c( "AIDS" ) ) {
      if( is.null( price ) ) {
         stop( "the 'AIDS' formula requires argument 'price'" )
      }
   }

   nGoods <- length( coef$alpha )
   nCoef  <- length( coef$all )

   createMatrix <- function( nGoods, nCoef, dim ) {
      result <- matrix( 0, nrow = nGoods^dim, ncol = nCoef )
      if( dim == 1 ) {
         rownames( result ) <- paste( "E", c( 1:nGoods ) )
      } else if( dim == 2 ) {
         rownames( result ) <- paste( "E", rep( 1:nGoods, each = nGoods ),
            rep( 1:nGoods, nGoods ) )
      }
      colnames( result ) <- names( coef$all )
      return( result )
   }

   jacobian <- list()
   jacobian$formula  <- formula
   jacobian$exp      <- createMatrix( nGoods, nCoef, 1 )
   jacobian$hicks    <- createMatrix( nGoods, nCoef, 2 )
   jacobian$marshall <- createMatrix( nGoods, nCoef, 2 )

   share <- array( share )

   aName <- paste( "alpha", c( 1:nGoods ) )
   bName <- paste( "beta", c( 1:nGoods ) )
   gName <- array( paste( "gamma", rep( 1:nGoods, nGoods ),
      rep( 1:nGoods, each = nGoods ) ), dim = c( nGoods, nGoods ) )
   eName <- array( paste( "E", rep( 1:nGoods, nGoods ),
      rep( 1:nGoods, each = nGoods ) ), dim = c( nGoods, nGoods ) )

   if( formula == "AIDS" ) {
      price <- array( price )
      for( i in 1:nGoods ) {
         # expenditure elasticities
         jacobian$exp[ paste( "E", i ), bName[ i ] ] <-
            1 / share[ i ]
         for( j in 1:nGoods ) {
            # Hicksian price elasticities
            jacobian$hicks[ eName[ i, j ], aName[ j ] ] <-
               -coef$beta[ i ] / share[ i ]
            jacobian$hicks[ eName[ i, j ], bName[ i ] ] <-
               - ( coef$alpha[ j ] - share[ j ] +
               coef$gamma[ j , ] %*% log( price ) ) / share[ i ]
            for( k in 1:nGoods ) {
               jacobian$hicks[ eName[ i, j ], gName[ k, j ] ] <-
                  ( i == k ) / share[ i ] -
                  coef$beta[ i ] * log( price[ k ] ) / share[ i ]
            }
            # Marshallian price elasticities
            jacobian$marshall[ eName[ i, j ], aName[ j ] ] <-
               -coef$beta[ i ] / share[ i ]
            jacobian$marshall[ eName[ i, j ], bName[ i ] ] <-
               - ( coef$alpha[ j ] +
               coef$gamma[ j , ] %*% log( price ) ) / share[ i ]
            for( k in 1:nGoods ) {
               jacobian$marshall[ eName[ i, j ], gName[ k, j ] ] <-
                  ( i == k ) / share[ i ] -
                  coef$beta[ i ] * log( price[ k ] ) / share[ i ]
            }
         }
      }
   } else {
      stop( "formula '", as.character( formula ), "' is not supported" )
   }
   return( jacobian )
}
