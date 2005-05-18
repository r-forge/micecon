invMillsRatio <- function( x, all = FALSE ) {
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
      library( VGAM )
      if( x@family@blurb[1] != "Bivariate probit model\n"  ) {
         stop( errorMessage )
      }
      library( mvtnorm )
      result <- data.frame( no = 1:nrow( x@predictors ),
         row.names = rownames( x@predictors ) )
      if( x@misc$link == "identity" ) {
         rho <- x@predictors[ , 3 ]
      } else if( x@misc$link == "rhobit" ){
         rho <- rhobit( x@predictors[ , 3 ], inv = TRUE )
      } else {
         stop( "the bivariate probit (binom2.rho) must be either estimated",
            " with link 'rhobit' or 'identity'" )
      }
      if( max( rho ) > 1 ) {
         stop( "the correlation between the error terms (rho) is larger",
            " than 1" )
      }
      pmvnormValues11 <- rep( NA, nrow( result ) )
      pmvnormValues10 <- rep( NA, nrow( result ) )
      pmvnormValues01 <- rep( NA, nrow( result ) )
      pmvnormValues00 <- rep( NA, nrow( result ) )
      for( i in 1:nrow( result ) ) {
         corrEq <- matrix( c( 1, rho[ i ], rho[ i ], 1 ), ncol = 2 )
         corrUneq <- matrix( c( 1, -rho[ i ], -rho[ i ], 1 ), ncol = 2 )
         if( x@y[ i, "11" ] == 1 | all ) {
            pmvnormValues11[ i ] <- pmvnorm( upper = x@predictors[ i, 1:2 ],
               corr = corrEq )
         }
         if( x@y[ i, "10" ] == 1 | all ) {
            pmvnormValues10[ i ] <- pmvnorm( upper = c( x@predictors[ i, 1 ],
               -x@predictors[ i, 2 ] ), corr = corrUneq )
         }
         if( x@y[ i, "01" ] == 1 | all ) {
            pmvnormValues01[ i ] <- pmvnorm( upper = c( -x@predictors[ i, 1 ],
               x@predictors[ i, 2 ] ), corr = corrUneq )
         }
         if( x@y[ i, "00" ] == 1 | all ) {
            pmvnormValues00[ i ] <- pmvnorm( upper = -x@predictors[ i, 1:2 ],
               corr = corrEq )
         }
      }

      result$IMR11a <- dnorm( x@predictors[ , 1 ] ) *
         pnorm( ( x@predictors[ , 2 ] - rho * x@predictors[ , 1 ] ) /
            ( 1 - rho^2 )^0.5 ) / pmvnormValues11
      result$IMR11b <- dnorm( x@predictors[ , 2 ] ) *
         pnorm( ( x@predictors[ , 1 ] - rho * x@predictors[ , 2 ] ) /
            ( 1 - rho^2 )^0.5 ) / pmvnormValues11
      result$IMR10a <- dnorm( x@predictors[ , 1 ] ) *
         pnorm( ( -x@predictors[ , 2 ] + rho * x@predictors[ , 1 ] ) /
            ( 1 - rho^2 )^0.5 ) / pmvnormValues10
      result$IMR10b <- -dnorm( x@predictors[ , 2 ] ) *
         pnorm( ( x@predictors[ , 1 ] - rho * x@predictors[ , 2 ] ) /
            ( 1 - rho^2 )^0.5 ) / pmvnormValues10
      result$IMR01a <- -dnorm( x@predictors[ , 1 ] ) *
         pnorm( ( x@predictors[ , 2 ] - rho * x@predictors[ , 1 ] ) /
            ( 1 - rho^2 )^0.5 ) / pmvnormValues01
      result$IMR01b <- dnorm( x@predictors[ , 2 ] ) *
         pnorm( ( -x@predictors[ , 1 ] + rho * x@predictors[ , 2 ] ) /
            ( 1 - rho^2 )^0.5 ) / pmvnormValues01
      result$IMR00a <- -dnorm( x@predictors[ , 1 ] ) *
         pnorm( ( -x@predictors[ , 2 ] + rho * x@predictors[ , 1 ] ) /
            ( 1 - rho^2 )^0.5 ) / pmvnormValues00
      result$IMR00b <- -dnorm( x@predictors[ , 2 ] ) *
         pnorm( ( -x@predictors[ , 1 ] + rho * x@predictors[ , 2 ] ) /
            ( 1 - rho^2 )^0.5 ) / pmvnormValues00
   } else {
      stop( errorMessage )
   }
   result$no <- NULL
   return( result )
}
