aidsConsist <- function( priceNames, totExpName, data, coef,
      priceIndex = "TL", basePrices = NULL, baseShares = NULL,
      shareNames = NULL ) {

   monoResult <- aidsMono( priceNames = priceNames, totExpName = totExpName,
      data = data, coef = coef, priceIndex = priceIndex,
      basePrices = basePrices, baseShares = baseShares )

   concResult <- aidsConcav( priceNames = priceNames, totExpName = totExpName,
      data = data, coef = coef, shareNames = shareNames )

   result <- list()
   result$mPercent <- monoResult$monoPercent
   result$monotony <- monoResult$monotony
   result$cPercent <- concResult$cPercent
   result$concavity <- concResult$concavity
   result$cMatrices <- concResult$cMatrices

   class( result ) <- "aidsConsist"
   return( result )
}
