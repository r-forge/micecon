checkConsist <- function( object, ... ) {
    UseMethod( "checkConsist" )
}

checkConsist.aidsEst <- function( object, ... ) {
   aidsConsist( priceNames = object$priceNames,
      shareNames = object$shareNames,
      totExpName = object$totExpName,
      data = get( as.character( object$call$data ) ),
      coef = object$coef,
      priceIndex = object$priceIndex,
      basePrices = object$basePrices,
      baseShares = object$baseShares )
}
