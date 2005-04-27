invMillsRatio <- function( x ) {
   errorMessage <- paste( "calculating the 'Inverse Mills Ratio' only works",
      "for probit models estimated by 'glm' or 'vglm'." )
   if( class( x )[ 1 ] == "glm" ) {
      if( x$family$family != "binomial" || x$family$link != "probit" ) {
         stop( errorMessage )
      }
      result <- data.frame( no = 1:nrow( x$data ),
         row.names = rownames( x$data ) )
      result$IMR1 <- dnorm( x$linear.predictors ) /
         pnorm( x$linear.predictors )
      result$delta1 <- result$IMR1 * ( result$IMR1 + x$linear.predictors )
      result$IMR0 <- dnorm( x$linear.predictors ) /
         ( 1 - pnorm( x$linear.predictors ) )
      result$delta0 <- result$IMR0 * ( result$IMR0 + x$linear.predictors )
   } else if( class( x ) == "vglm" ) {
      if( x@family@blurb[1] != "Bivariate probit model\n"  ) {
         stop( errorMessage )
      }
      library( mvtnorm )
      result <- data.frame( no = 1:nrow( x@predictors ),
         row.names = rownames( x@predictors ) )
      ya <- rowSums( x@y[ , c( "10", "11" ) ] )
      yb <- rowSums( x@y[ , c( "01", "11" ) ] )
      rho <- x@predictors[ , 3 ]
      pmvnormValues <- rep( NA, nrow( result ) )
      for( i in seq( along = pmvnormValues ) ) {
         corr <- matrix( c( 1, rho[ i ], rho[ i ], 1 ), ncol = 2 )
         pmvnormValues[ i ] <- pmvnorm( upper = - x@predictors[ i, 1:2 ],
            corr = corr )
      }
      result$IMRa1 <- dnorm( x@predictors[ , 1 ] ) *
         pnorm( ( - x@predictors[ , 2 ] - rho * ya ) / ( 1 - rho^2 )^0.5 ) /
         pmvnormValues
      result$IMR1 <- dnorm( x@predictors[ , 1 ] ) /
         pnorm( x@predictors[ , 1 ] )
   } else {
      stop( errorMessage )
   }
   return( result )
}

