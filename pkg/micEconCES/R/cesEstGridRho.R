cesEstGridRho <- function( from = -0.8, to = 4, by = 0.2,
      returnAll = FALSE, ... )  {

   # some tests
   if( from < -1 ) {
      stop( "argument 'from' must not be smaller than '-1'" )
   } else if( from > by ) {
      stop( "argument 'from' must be larger than argument 'to'" )
   } else if( by <= 0 ) {
      stop( "argument 'by' must be positive" )
   }

   # list that should contain each estimation result
   allResults <- list()

   # summary results for each estimation (with different fixed rhos)
   sumResults <- data.frame( rho = seq( from = from, to = to, by = by ) )
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
   result$otherSum <- sumResults

   # add full results of each estimation
   if( returnAll ) {
      result$otherFull <- allResults
   }

   class( result ) <- c( "cesEstRhoGrid", class( result ) )
   return( result )
}
