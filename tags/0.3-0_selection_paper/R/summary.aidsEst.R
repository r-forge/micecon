summary.aidsEst <- function( object, elaFormula = NULL,
      quantNames = NULL, ... ) {

   result <- object

   # specify default value for argument elaFormula
   if( is.null ( elaFormula ) ) {
      if( substr( object$method, 1, 2 ) == "LA" ) {
         elaFormula <- "Ch"
      } else if( substr( object$method, 1, 2 ) %in% c( "MK", "IL" ) ) {
         elaFormula <- "AIDS"
      }
   }

   # test reasonability of argument elaFormula
   if( substr( object$method, 1, 2 ) %in% c( "MK", "IL" ) &&
         elaFormula != "AIDS" ) {
      warning( paste( "It does not make sense to calculate the elasticities",
         " of a (non-linear) AIDS model with formula '", elaFormula, "'",
         sep = "" ) )
   }

   # test argument quantNames
   if( !is.null( quantNames ) &&
         length( quantNames ) != length( object$wMeans ) ){
      stop( "argument 'quantNames' must be a vector with as many elements",
         " as there are goods (in this case ", length( object$wMeans ), ")" )
   }

   # to avoid warning message in aidsEla
   if( elaFormula %in% c( "Ch", "EU" ) ) {
      object$pMeans <- NULL
   }

   # calculate demand elasticities
   result$ela  <- aidsEla( coef = result$coef,
      shares = object$wMeans, prices = object$pMeans,
      formula = elaFormula,
      priceNames = object$priceNames, quantNames = quantNames,
      coefVcov = object$coef$allcov, df = object$est$df )

   result$elaFormula <- elaFormula

   class( result ) <- "summary.aidsEst"
   return( result )
}

