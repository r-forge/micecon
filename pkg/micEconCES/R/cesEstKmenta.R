cesEstKmenta <- function( yName, xNames, data ){

   result <- list()

   ## Estimating unrestricted model
   result$translogUnr <- translogEst( yName = yName, xNames = xNames, data = data )

   ## Testing linear approximation - Wald test
   result$test <- lht( result$translogUnr$est,
      c( "b_1_2 = -b_1_1", "b_1_2 = -b_2_2" ))

   ## Estimating restricted model
   result$translog <- systemfit( formula = formula(result$translog$est),
      data = model.frame( result$translog$est ),
      restrict.matrix = c( "eq1_b_1_2 = -eq1_b_1_1", "eq1_b_1_2 = -eq1_b_2_2" ))

   ## Parameter vector
   result$coefficients <- numeric( 4 )
   names( result$coefficients ) <- c( "gamma", "delta", "rho", "phi" )

   ## Defining gamma
   result$coefficients[ "gamma" ] <- exp( coef( result$translog )[ "eq1_(Intercept)" ] )

   ## Defining phi
   result$coefficients[ "phi" ] <- coef( result$translog )[ "eq1_a_1" ] +
                             coef( result$translog )[ "eq1_a_2" ]

   ## Defining delta
   result$coefficients[ "delta" ] <- coef( result$translog )[ "eq1_a_1" ] /
                               result$coefficients[ "phi" ]

   ## Defining rho
   result$coefficients[ "rho" ] <- coef( result$translog )[ "eq1_b_1_2" ] /
      ( coef( result$translog )[ "eq1_a_1" ] * coef( result$translog )[ "eq1_a_2" ] /
         result$coefficients[ "phi" ] )

   ## Delta method
   jacobian <- matrix( 0, nrow = length( result$coefficients ),
      ncol = length( coef( result$translog )))
   rownames( jacobian ) <- names( result$coefficients )
   colnames( jacobian ) <- names( coef( result$translog ))

   jacobian[ "gamma", "eq1_(Intercept)" ] <-
      exp( coef( result$translog )[ "eq1_(Intercept)" ] )
   jacobian[ "phi", c( "eq1_a_1", "eq1_a_2" )] <- 1
   jacobian[ "delta", "eq1_a_1" ] <- 1 / result$coefficients[ "phi" ] -
      coef( result$translog )[ "eq1_a_1" ] / result$coefficients[ "phi" ]^2
   jacobian[ "delta", "eq1_a_2" ] <-
      - coef( result$translog )[ "eq1_a_1" ] / result$coefficients[ "phi" ]^2
   jacobian[ "rho", "eq1_a_1" ] <- coef( result$translog )[ "eq1_b_1_2" ] /
      coef( result$translog )[ "eq1_a_1" ] * coef( result$translog )[ "eq1_a_2" ] -
      coef( result$translog )[ "eq1_b_1_2" ] * result$coefficients[ "phi" ] /
      coef( result$translog )[ "eq1_a_1" ]^2 / coef( result$translog )[ "eq1_a_2" ]^3
   jacobian[ "rho", "eq1_a_2" ] <- coef( result$translog )[ "eq1_b_1_2" ] /
      coef( result$translog )[ "eq1_a_1" ] * coef( result$translog )[ "eq1_a_2" ] -
      coef( result$translog )[ "eq1_b_1_2" ] * result$coefficients[ "phi" ] /
      coef( result$translog )[ "eq1_a_1" ]^3 / coef( result$translog )[ "eq1_a_2" ]^2
   jacobian[ "rho", "eq1_b_1_2" ] <-  result$coefficients[ "phi" ] /
      coef( result$translog )[ "eq1_a_1" ] / coef( result$translog )[ "eq1_a_2" ]
   result$vcov <- jacobian %*% vcov( result$translog ) %*%  t( jacobian )

   return(result)
}

