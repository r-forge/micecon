quadFuncEst <- function( yName, xNames, data, shifterNames = NULL,
   quadHalf = TRUE, regScale = 1, ... ) {

   checkNames( c( yName, xNames, shifterNames ), names( data ) )

   nExog   <- length( xNames )
   nShifter <- length( shifterNames )
   result <- list()

   if( "plm.dim" %in% class( data ) ) {
      estData <- data[ , 1:2 ]
      estData$y <- data[[ yName ]]
   } else {
      estData <- data.frame( y = data[[  yName ]] )
   }

   estFormula <- "y ~ 1"
   if( nExog > 0 ) {
      for( i in 1:nExog ) {
         xName <- paste( "x", as.character( i ), sep = "" )
         estData[[ xName ]] <- data[[ xNames[ i ] ]] / regScale
         estFormula <- paste( estFormula, "+", xName )
      }
      for( i in 1:nExog ) {
         for( j in i:nExog ) {
            xName <- paste( "x", as.character( i ), ".", as.character( j ),
               sep = "" )
            estData[[ xName ]] <- ifelse( quadHalf, 0.5, 1 ) *
               ifelse( i == j, 1, 2 ) *
               data[[ xNames[ i ] ]] * data[[ xNames[ j ] ]] / regScale
            estFormula <- paste( estFormula, "+", xName )
         }
      }
   }
   if( nShifter > 0 ) {
      for( i in 1:nShifter ) {
         if( is.factor( data[[ shifterNames[ i ] ]] ) | 
               is.logical( data[[ shifterNames[ i ] ]] ) ) {
            xName <- paste( "s", as.character( i ), ".", sep = "" )
            estData[[ xName ]] <- data[[ shifterNames[ i ] ]]
         } else {
            xName <- paste( "s", as.character( i ), sep = "" )
            estData[[ xName ]] <- data[[ shifterNames[ i ] ]] / regScale
         }
         estFormula <- paste( estFormula, "+", xName )
      }
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
   coefNames <- .quadFuncCoefNames( nExog, nShifter )
   for( i in names( result$est$xlevels ) ) {
      shiftNum <- as.integer( sub( "^s", "", sub( "\\.$", "", i ) ) )
      coefNum <- which( coefNames == paste( "d", shiftNum, sep = "_" ) )
      coefNamesFac <- grep( paste( "^", i, sep = "" ), names( result$coef ), 
         value = TRUE )
      newCoefNames <- paste( "d", shiftNum, 
         sub( paste( "^", i, sep = "" ), "", coefNamesFac ),
         sep = "_" )
      coefNames <- c( coefNames[ 1:( coefNum - 1 ) ], newCoefNames,
         coefNames[ -c( 1:coefNum ) ] )
   }
   names( result$coef )      <- coefNames
   rownames( result$coefCov ) <- coefNames
   colnames( result$coefCov ) <- coefNames

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
