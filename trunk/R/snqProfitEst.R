
## ========= Estimation =============
snqProfitEst <- function( pNames, qNames, fNames = NULL,
   ivNames = NULL, data,  form = 0, base = 1,
   weights = snqProfitWeights( pNames, qNames, data, "DW92", base = base ),
   method = ifelse( is.null( ivNames ), "SUR", "3SLS" ), ... ) {

   if( length( qNames ) != length( pNames ) ) {
      stop( "arguments 'qNames' and 'pNames' must have the same length." )
   }
   if( length( pNames ) < 2 ) {
      stop( "you must specify at least 2 netputs." )
   }
   if( length( pNames ) != length( weights ) ) {
      stop( "arguments 'pNames' and 'weights' must have the same length." )
   }
   if( min( weights ) < 0 ) {
      warning( paste( "At least one weight of the prices for normalization",
         "(argument 'weights') is negative. Thus, in this case positive",
         "semidefiniteness of the 'beta' matrix does not ensure",
         "a convex profit function." ) )
   }

   nNetput <- length( qNames )  # number of netputs
   nFix    <- length( fNames )  # number of fixed inputs
   nIV     <- length( ivNames )  # number of fixed inputs
   nObs    <- nrow( data )      # number of observations

   if( form == 0 ) {
      nCoef   <- nNetput + nNetput * ( nNetput - 1 ) / 2 + nNetput * nFix +
         ( nFix + 1 ) * nFix/2  #number of coefficients
   } else if( form == 1 ) {
      nCoef   <- nNetput + nNetput * ( nNetput - 1 ) / 2 + nNetput * nFix +
         nNetput * ( nFix + 1 ) * nFix/2  #number of coefficients
   } else {
      stop( "argument 'form' must be either 0 or 1." )
   }

   result  <- list()

   ## scaling data
   if( !is.null( base ) ) {
      scaledData <- data.frame( nr = 1:nObs )
      for( i in 1:nNetput ) {
         scaledData[[ pNames[ i ] ]] <- data[[ pNames[ i ] ]] /
            mean( data[[ pNames[ i ] ]][ base ] )
         scaledData[[ qNames[ i ] ]] <- data[[ qNames[ i ] ]] *
            mean( data[[ pNames[ i ] ]][ base ] )
      }
   } else {
      scaledData <- data
   }

   ## price index for normalization
   estData <- data.frame( nr = 1:nObs, normPrice = 0 )
   for( i in 1:nNetput ) {
      estData$normPrice <- estData$normPrice +
         with( scaledData, get( pNames[ i ] ) ) * weights[ i ]
   }

   ## real/normalized netput prices and netput quantities
   result$pMeans <- array( NA, nNetput )
   result$qMeans <- array( NA, nNetput )
   for( i in 1:nNetput ) {
      estData[[ paste( "pr", as.character( i ), sep = "" ) ]] <-
         with( scaledData, get( pNames[ i ] ) ) / estData$normPrice
      result$pMeans[ i ] <- mean( with( scaledData, get( pNames[ i ] ) ) )
      estData[[ paste( "q", as.character( i ), sep = "" ) ]] <-
         with( scaledData, get( qNames[ i ] ) )
      result$qMeans[ i ] <- mean( with( scaledData, get( qNames[ i ] ) ) )
   }
   names( result$pMeans ) <- pNames
   names( result$qMeans ) <- qNames

   ## quadratic netput prices
   for( i in 1:nNetput ) {
      for( j in 1:nNetput ) {
         for( k in 1:nNetput ) {
            estData[[ paste( "pq", as.character( i ), ".", as.character( j ), ".",
               as.character( k ), sep = "" ) ]] <-
               -0.5 * weights[ i ] * with( scaledData, get( pNames[ j ] ) ) *
               with( scaledData, get( pNames[ k ] ) ) / estData$normPrice^2
         }
      }
   }

   ## quasi-fix inputs
   if( nFix > 0 ) {
      for( i in 1:nFix ) {
         estData[[ paste( "f", as.character( i ), sep = "" ) ]] <-
            with( data, get( fNames[ i ] ) )
         result$fMeans[ i ] <- mean( with( data, get( fNames[ i ] ) ) )
      }
      ## quadratic quasi-fix inputs
      for( i in 1:nNetput ) {
         for( j in 1:nFix ) {
            for( k in 1:nFix ) {
               estData[[ paste( "fq", as.character( i ), ".", as.character( j ), ".",
               as.character( k ), sep = "" ) ]] <-
                  0.5 * weights[ i ] * with( data, get( fNames[ j ] ) ) *
                  with( data, get( fNames[ k ] ) )
            }
         }
      }
   }

   ## instrumental variables
   if( nIV == 0 ) {
      inst <- NULL
   } else {
      inst <- as.formula( paste( "~", paste( ivNames, collapse = "+" ) ) )
      for( i in 1:nIV ) {
         estData[[ paste( "iv", as.character( i ), sep = "" ) ]] <-
            with( data, get( ivNames[ i ] ) )
      }
   }
   system <- snqProfitSystem( nNetput, nFix )    # equation system
   restrict <- snqProfitRestrict( nNetput, nFix, form )    # restrictions
   result$est <- systemfit( method = method, eqns = system, data = estData,
      TX = restrict, inst = inst, ... )
   result$coef <- snqProfitCoef( coef = result$est$bt, nNetput = nNetput,
      nFix = nFix, form = form, coefCov = result$est$btcov,
      df = nNetput * nObs - nCoef )
      # estimated coefficients
   result$r2 <- array( NA, c( nNetput ) )
   for( i in 1:nNetput ) {
      result$r2[ i ] <- result$est$eq[[ i ]]$r2
   }

   result$hessian <- snqProfitHessian( result$coef$beta, result$pMeans, weights )
      # Hessian matrix
   result$ela <- snqProfitEla( result$coef$beta, result$pMeans,
      result$qMeans, weights )   # estimated elasticities
   result$scaledData  <- scaledData
   result$coef$liCoef <- result$est$bt
   result$coef$liCoefCov <- result$est$btcov
   result$weights  <- weights
   result$normPrice <- estData$normPrice
   result$convexity  <- semidefiniteness( result$hessian[
      1:( nNetput - 1 ), 1:( nNetput - 1 ) ] )$positive
   class( result )  <- "snqProfitEst"
   return( result )
}
