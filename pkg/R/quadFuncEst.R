quadFuncEst <- function( yName, xNames, data, shifterNames = NULL,
   linear = FALSE, homWeights = NULL, quadHalf = TRUE, regScale = 1, ... ) {

   checkNames( c( yName, xNames, shifterNames ), names( data ) )

   # check argument 'homWeights'
   if( !is.null( homWeights ) ) {
      if( is.null( names( homWeights ) ) ) {
         stop( "the elements of argument 'homWeights' must have names" )
      }
      if( !all( names( homWeights ) %in% xNames ) ) {
         stop( "all names in argument 'homWeights' must be in argument 'xNames'" )
      }
      if( abs( sum( homWeights ) - 1 ) > .Machine$double.eps ^ 0.5 ) {
         stop( "the sum of the elements in argument 'homWeights' must be 1" )
      }
   }

   nExog   <- length( xNames )
   nShifter <- length( shifterNames )
   result <- list()
   result$call <- match.call()

   if( "plm.dim" %in% class( data ) ) {
      estData <- data[ , 1:2 ]
      estData$y <- data[[ yName ]]
   } else {
      estData <- data.frame( y = data[[  yName ]] )
   }

   if( !is.null( homWeights ) ) {
      if( ! linear ) {
         stop( "imposing homogeneity of degree zero currently works",
            "only for linear functions" )
      }
      estData$deflator <- 0
      for( i in seq( along = homWeights ) ) {
         estData$deflator <- estData$deflator + 
            homWeights[ i ] * data[[ names( homWeights )[ i ] ]]
      }
      whichHom <- which( xNames %in% names( homWeights ) )
      iOmit <- which( xNames == names( homWeights )[ length( homWeights ) ] )
   } else {
      iOmit <- 0
   }

   estFormula <- "y ~ 1"
   for( i in seq( along = xNames ) ) {
      if( i != iOmit ) {
         xName <- paste( "a", as.character( i ), sep = "_" )
         estData[[ xName ]] <- .quadFuncVarHom( data, xNames[ i ], 
            homWeights, estData$deflator ) / regScale
         estFormula <- paste( estFormula, "+", xName )
      }
   }
   if( !linear ) {
      for( i in seq( along = xNames ) ) {
         for( j in i:nExog ) {
            xName <- paste( "b", as.character( i ), as.character( j ),
               sep = "_" )
            estData[[ xName ]] <- ifelse( quadHalf, 0.5, 1 ) *
               ifelse( i == j, 1, 2 ) *
               data[[ xNames[ i ] ]] * data[[ xNames[ j ] ]] / regScale
            estFormula <- paste( estFormula, "+", xName )
         }
      }
   }
   for( i in seq( along = shifterNames ) ) {
      if( is.factor( data[[ shifterNames[ i ] ]] ) | 
            is.logical( data[[ shifterNames[ i ] ]] ) ) {
         xName <- paste( "d", "_", as.character( i ), "_", sep = "" )
         estData[[ xName ]] <- data[[ shifterNames[ i ] ]]
      } else {
         xName <- paste( "d", as.character( i ), sep = "_" )
         estData[[ xName ]] <- data[[ shifterNames[ i ] ]] / regScale
      }
      estFormula <- paste( estFormula, "+", xName )
   }
   result$nExog <- nExog
   result$nShifter <- nShifter
   if( "plm.dim" %in% class( data ) ) {
      result$est <- plm( as.formula( estFormula ), estData, ... )
      result$est$call$formula <- as.formula( estFormula )
   } else {
      result$est <- lm( as.formula( estFormula ), estData, ... )
   }
   result$residuals <- residuals( result$est )
   result$fitted    <- estData$y - result$residuals

   # coefficients and their covariance matrix
   result$coef      <- coef( result$est )
   result$coefCov   <- vcov( result$est )
   if( "plm.dim" %in% class( data ) ) {
      if( is.null( result$est$call$model ) ||
            result$est$call$model == "within" ) {
         result$coef <- c( result$est$alpha, result$coef )
         result$coefCov <- rbind( NA, cbind( NA, vcov( result$est ) ) )
      }
   }
   names( result$coef )[ 1 ]       <- "a_0"
   rownames( result$coefCov )[ 1 ] <- "a_0"
   colnames( result$coefCov )[ 1 ] <- "a_0"

   # adding coefficient that has been dropped due to the homogeneity restriction
   # and its covariances
   if( !is.null( homWeights ) ) {
      # missing coefficient
      coefOmit <- 0
      for( i in whichHom[ whichHom != iOmit ] ) {
         coefOmit <- coefOmit - result$coef[ paste( "a", i, sep = "_" ) ]
      }
      result$coef <- c( result$coef[ 1:iOmit ], coefOmit, 
         result$coef[ -c( 1:iOmit ) ] )
      names( result$coef )[ iOmit + 1 ] <- paste( "a", iOmit, sep = "_" )
      # missing row of covariance matrix
      coefCovOmit <- rep( 0, ncol( result$coefCov ) )
      for( i in whichHom[ whichHom != iOmit ] ) {
         coefCovOmit <- coefCovOmit - 
            result$coefCov[ paste( "a", i, sep = "_" ), ]
      }
      result$coefCov <- rbind( result$coefCov[ 1:iOmit, ], coefCovOmit, 
         result$coefCov[ -c( 1:iOmit ), ] )
      rownames( result$coefCov )[ iOmit + 1 ] <- paste( "a", iOmit, sep = "_" )
      # missing column of covariance matrix
      coefCovOmit <- rep( 0, nrow( result$coefCov ) )
      for( i in whichHom[ whichHom != iOmit ] ) {
         coefCovOmit <- coefCovOmit - 
            result$coefCov[ , paste( "a", i, sep = "_" ) ]
      }
      result$coefCov <- cbind( result$coefCov[ , 1:iOmit ], coefCovOmit, 
         result$coefCov[ , -c( 1:iOmit ) ] )
      colnames( result$coefCov )[ iOmit + 1 ] <- paste( "a", iOmit, sep = "_" )
   }

   if( linear & nExog > 0 ) {
      nQuadCoef <- nExog * ( nExog + 1 ) / 2
      quadCoefNames <- paste( "b", 
         vecli( matrix( rep( 1:nExog, nExog ), nrow = nExog ) ), 
         vecli( matrix( rep( 1:nExog, each = nExog ), nrow = nExog ) ),
         sep = "_" )
      quadCoef <- rep( 0, nQuadCoef )
      names( quadCoef ) <- quadCoefNames
      result$coef <- c( result$coef[ 1:( nExog + 1 ) ], quadCoef,
         result$coef[ -c( 1:( nExog + 1 ) ) ] )
      quadCoefCovRows <- matrix( 0, nrow = nQuadCoef, 
         ncol = ncol( result$coefCov ) ) 
      rownames( quadCoefCovRows ) <- quadCoefNames
      result$coefCov <- rbind( result$coefCov[ 1:( nExog + 1 ), ],
         quadCoefCovRows, result$coefCov[ -c( 1:( nExog + 1 ) ), ] )
      quadCoefCovCols <- matrix( 0, nrow = nrow( result$coefCov ), 
         ncol = nQuadCoef ) 
      colnames( quadCoefCovCols ) <- quadCoefNames
      result$coefCov <- cbind( result$coefCov[ , 1:( nExog + 1 ) ],
         quadCoefCovCols, result$coefCov[ , -c( 1:( nExog + 1 ) ) ] )
   }

   result$r2    <- summary( result$est )$r.squared
   result$r2bar <- summary( result$est )$adj.r.squared
   result$nObs  <- length( result$residuals )
   if( "plm.dim" %in% class( data ) ) {
      result$model.matrix <- cbind( rep( 1, result$nObs ),
         as.matrix( estData[ , 4:( ncol( estData ) ) ] ) )
   } else {
      result$model.matrix <- cbind( rep( 1, result$nObs ),
         as.matrix( estData[ , 2:( ncol( estData ) ) ] ) )
   }
   class( result ) <- "quadFuncEst"
   return( result )
}
