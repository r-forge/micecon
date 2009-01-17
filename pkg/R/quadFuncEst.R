quadFuncEst <- function( yName, xNames, data, shifterNames = NULL,
   quadHalf = TRUE, regScale = 1, ... ) {

   checkNames( c( yName, xNames, shifterNames ), names( data ) )

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

   estFormula <- "y ~ 1"
   if( nExog > 0 ) {
      for( i in 1:nExog ) {
         xName <- paste( "a", as.character( i ), sep = "_" )
         estData[[ xName ]] <- data[[ xNames[ i ] ]] / regScale
         estFormula <- paste( estFormula, "+", xName )
      }
      for( i in 1:nExog ) {
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
   if( nShifter > 0 ) {
      for( i in 1:nShifter ) {
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
