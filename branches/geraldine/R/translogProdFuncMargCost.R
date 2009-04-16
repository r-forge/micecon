translogProdFuncMargCost <- function( yName, xNames, wNames,
      data, coef, dataLogged = FALSE ) {

   checkNames( c( yName, xNames, wNames ), names( data ) )

   if( length( yName ) != 1 ) {
      stop( "argument 'yName' must include only one input" )
   }

   if( length( xNames ) != length( wNames ) ) {
      stop( "arguments 'xNames' and 'wNames' must have the same length" )
   }

   if( dataLogged ) {
      logData   <- data
   } else {
      logData <- .micEconLogData( data = data,
         varNames = c( yName, xNames, wNames ) )
   }

   # number of inputs
   nInput <- length( xNames )

   # elasticities: d log F / d log x_i
   elaData <- translogEla( xNames = xNames, data = logData, coef = coef,
      dataLogged = TRUE )

   # matrix of beta coefficients
   beta <- matrix( NA, nrow = nInput, ncol = nInput )
   for( i in 1:nInput ) {
      for( j in i:nInput ) {
         beta[ i, j ] <- coef[ paste( "b", i, j, sep = "_" ) ]
         beta[ j, i ] <- beta[ i, j ]
      }
   }

   # marginal costs
   margCost <- matrix( NA, nrow = nrow( logData ), ncol = 1 )

   for( t in 1:nrow( logData ) ) {
      yVal <- exp( as.numeric( logData[ t, yName ] ) )
      xVal <- exp( as.numeric( logData[ t, xNames ] ) )
      wVal <- exp( as.numeric( logData[ t, wNames ] ) )
      eVal <- as.numeric( elaData[ t, ] )

      # Jacobian matrix of g with respect to x
      gxJac <- matrix( NA, nrow = nInput, ncol = nInput )
      for( j in 1:nInput ) {
         for( i in 1:( nInput - 1 ) ) {
            gxJac[ i, j ] <-
               ( i == j ) * wVal[ nInput ] *
                  ( xVal[ nInput ] / xVal[ j ]^2 ) *
                  eVal[ i ] / eVal[ nInput ] -
               ( j == nInput ) * wVal[ nInput ] * ( 1 / xVal[ i ] ) *
                  eVal[ i ] / eVal[ nInput ] -
               wVal[ nInput ] * ( xVal[ nInput ] / xVal[ i ] ) *
                  ( beta[ i, j ] / xVal[ j ] ) / eVal[ nInput ] +
               wVal[ nInput ] * ( xVal[ nInput ] / xVal[ i ] ) *
                  ( eVal[ i ] / eVal[ nInput ]^2 ) *
                  beta[ nInput, j ] / xVal[ j ]
         }
         gxJac[ nInput, j ] <- eVal[ j ] / xVal[ j ]
      }

      # Jacobian matrix of g with respect to y
      gyJac <- matrix( nrow = nInput, ncol = 1 )
      gyJac[ 1:( nInput - 1 ), 1 ] <- 0
      gyJac[ nInput, 1 ] <- - 1 / yVal

      # Jacobian matrix of x with respect to y
      xyJac <- - solve( gxJac, gyJac )

      # marginal costs
      margCost[ t, ] <- t( wVal ) %*% xyJac
   }

   margCost <- as.data.frame( margCost )
   rownames( margCost ) <- rownames( data )
   names( margCost ) <- yName

   return( margCost )
}