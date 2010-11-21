cesEstGridRho <- function( rho1Values, rhoValues, returnAll, ... )  {

   # some tests
   if( !is.null( rho1Values ) && !is.null( rhoValues ) ) {
      stop( "currently, either argument 'rho1Values' or argument 'rhoValues'",
         " must be NULL" )
   }
   if( is.null( rho1Values ) && is.null( rhoValues ) ) {
      stop( "either argument 'rho1Values' or argument 'rhoValues'",
         " must be non-NULL" )
   }
   if( !is.null( rho1Values ) ) {
      if( !is.numeric( rho1Values ) ) {
         stop( "the rho_1s specified in argument 'rho1Values'",
            " must be numeric" )
      } else if(  min( rho1Values ) < -1 ) {
         stop( "the rho_1s specified in argument 'rho1Values'",
            " must not be smaller than '-1'" )
      }
   }
   if( !is.null( rhoValues ) ) {
      if( !is.numeric( rhoValues ) ) {
         stop( "the rhos specified in argument 'rhoValues'",
            " must be numeric" )
      } else if( min( rhoValues ) < -1 ) {
         stop( "the rhos specified in argument 'rhoValues'",
            " must not be smaller than '-1'" )
      }
   }

   # list that should contain each estimation result
   allResults <- list()

   # summary results for each estimation (with different fixed rhos)
   if( !is.null( rho1Values ) ) {
      sumResults <- data.frame( rho1 = rho1Values )
   } else {
      sumResults <- data.frame( rho = rhoValues )
   }
   sumResults$rss <- NA
   sumResults$convergence <- NA

   # estimate the CES for each pre-defined rho
   for( i in 1:nrow( sumResults ) ) {
      allResults[[ i ]] <- cesEst( rho1 = sumResults[[ "rho1" ]][ i ], 
         rho = sumResults[[ "rho" ]][ i ], ... )
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
