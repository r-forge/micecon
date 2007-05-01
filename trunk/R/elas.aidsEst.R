elas.aidsEst <- function( object, formula = NULL, ... ) {

   # specify default value for argument formula
   if( is.null ( formula ) ) {
      if( substr( object$method, 1, 2 ) == "LA" ) {
         formula <- "Ch"
      } else if( substr( object$method, 1, 2 ) %in% c( "MK", "IL" ) ) {
         formula <- "AIDS"
      }
   }

   # test reasonability of argument formula
   if( substr( object$method, 1, 2 ) %in% c( "MK", "IL" ) &&
         formula != "AIDS" ) {
      warning( paste( "It does not make sense to calculate the elasticities",
         " of a (non-linear) AIDS model with formula '", formula, "'",
         sep = "" ) )
   }

   # to avoid warning message in aidsElas
   if( formula %in% c( "Ch", "EU" ) ) {
      object$pMeans <- NULL
   }

   # calculate demand elasticities
   result  <- aidsElas( coef = object$coef,
      shares = object$wMeans, prices = object$pMeans,
      formula = formula,
      priceNames = object$priceNames,
      coefVcov = object$coef$allcov, df = object$est$df, ... )

   return( result )
}

