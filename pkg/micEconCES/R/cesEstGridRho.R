cesEstGridRho <- function( rhoValues, returnAll, ... )  {

   # some tests
   if( min( rhoValues ) < -1 ) {
      stop( "the rhos specified in argument 'rhoValues'",
         " must not be smaller than '-1'" )
   }

   # list that should contain each estimation result
   allResults <- list()

   # summary results for each estimation (with different fixed rhos)
   sumResults <- data.frame( rho = rhoValues )
   sumResults$rss <- NA
   sumResults$convergence <- NA

   # estimate the CES for each pre-defined rho
   for( i in 1:nrow( sumResults ) ) {
      allResults[[ i ]] <- cesEst( rho = sumResults$rho[ i ], ... )
      sumResults$rss[ i ] <- allResults[[ i ]]$rss
      if( !is.null( allResults[[ i ]]$convergence ) ) {
         sumResults$convergence[ i ] <- allResults[[ i ]]$convergence
      }
   }

   # returned object: the estimation results with the lowest RSS
   result <- allResults[[ which.min( sumResults$rss ) ]]

   # add the summary results of each estimation
   result$allRhoSum <- sumResults

   # add full results of each estimation
   if( returnAll ) {
      result$allRhoFull <- allResults
   }

   result$call <- match.call()
   return( result )
}
