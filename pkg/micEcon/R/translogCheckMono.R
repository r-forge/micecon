translogCheckMono <- function( xNames, data, coef, increasing = TRUE,
   strict = FALSE, dataLogged = FALSE,
   tol = 10 * .Machine$double.eps ) {

   result <- list()

   deriv <- translogDeriv( xNames = xNames, data = data, coef = coef,
      dataLogged = dataLogged )$deriv

   nExog <- ncol( deriv )
   nObs <- nrow( deriv )

   result$exog <- matrix( NA, nrow = nObs, ncol = nExog )
   colnames( result$exog ) <- colnames( deriv )

   for( i in 1:nExog ) {
      if( strict ) {
         result$exog[ , i ] <- deriv[ , i ] * (-1)^increasing < 0
      } else {
         result$exog[ , i ] <- deriv[ , i ] * (-1)^increasing <= tol
      }
   }

   result$obs <- rowSums( !result$exog ) == 0
   result$increasing <- increasing
   result$strict     <- strict

   class( result ) <- "translogCheckMono"
   return( result )
}

