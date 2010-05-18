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
      if( increasing ) {
         if( strict ) {
            result$exog[ , i ] <- deriv[ , i ] > 0
         } else {
            result$exog[ , i ] <- deriv[ , i ] >= - tol
         }
      } else {
         if( strict ) {
            result$exog[ , i ] <- deriv[ , i ] < 0
         } else {
            result$exog[ , i ] <- deriv[ , i ] <= tol
         }
      }
   }

   result$obs <- rowSums( !result$exog ) == 0
   result$increasing <- increasing
   result$strict     <- strict

   class( result ) <- "translogCheckMono"
   return( result )
}

