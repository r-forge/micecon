quadFuncDeriv <- function( xNames, data, coef, coefCov = NULL,
   quadHalf = TRUE ) {

   checkNames( c( xNames ), names( data ) )

   result <- list()

   nExog <- length( xNames )
   nCoef <- 1 + nExog + nExog * ( nExog + 1 ) / 2

   if( nCoef > length( coef ) ) {
      stop( "a quadratic function with ", nExog, " exogenous variables",
         " must have at least ", nCoef, " coefficients" )
   }

   ## derivatives
   deriv <- array( NA, c( nrow( data ), nExog ) )
   for( i in 1:nExog ) {
      deriv[ , i ] <- coef[ paste( "a", i, sep = "_" ) ]
      for( j in 1:nExog ) {
         deriv[ , i ] <- deriv[ , i ] + ifelse( quadHalf, 1, 2 ) *
            coef[ paste( "b", min( i, j ), max( i, j ), sep = "_" ) ] * 
            data[[ xNames[ j ] ]]
      }
   }
   colnames( deriv ) <- xNames
   result$deriv    <- as.data.frame( deriv )

   if( !is.null( coefCov ) ) {
      ## variances of the derivatives
      variance <- array( NA, c( nrow( data ), nExog ) )
      for(i in 1:nExog ) {
         variance[ , i ] <- coefCov[ paste( "a", i, sep = "_" ), 
            paste( "a", i, sep = "_" ) ]   # variance of aplha(i)
         for( j in 1:nExog ) {
            variance[ , i ] <- variance[ , i ] +
               coefCov[ paste( "a", i, sep = "_" ), 
                  paste( "b", min( i, j ), max( i, j ), sep = "_" ) ] *
               ifelse( quadHalf, 1, 2 ) * data[[ xNames[ j ] ]]
               # covariance alpha(i)-beta(i,_)
         }
         for( j in 1:nExog ) {
            for( k in 1:nExog ) {
               variance[ , i ] <- variance[ , i ] +
                  coefCov[ paste( "b", min( i, j ), max( i, j ), sep = "_" ),
                     paste( "b", min( i, k ), max( i, k ), sep = "_" ) ] *
                  ifelse( quadHalf, 1, 4 ) *
                  data[[ xNames[ j ] ]] * data[[ xNames[ k ] ]]
                  # variances + covariance beta(i,_)-beta(i,_)
            }
         }
      }
      stdDev <- variance^0.5  # standard errors
      colnames( variance ) <- xNames
      colnames( stdDev )   <- xNames
      result$variance <- as.data.frame( variance )
      result$stdDev   <- as.data.frame( stdDev )
   }

   class( result ) <- "quadFuncDeriv"
   return( result )
}
