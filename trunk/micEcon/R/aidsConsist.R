aidsConsist <- function( priceNames, totExpName, data, coef,
      priceIndex = "TL", basePrices = NULL, baseShares = NULL,
      shareNames = NULL ) {

   if( priceIndex == "TL" ){
      shareNamesMono <- NULL
   } else {
      shareNamesMono <- shareNames
   }

   monoResult <- aidsMono( priceNames = priceNames, totExpName = totExpName,
      data = data, coef = coef, priceIndex = priceIndex,
      basePrices = basePrices, baseShares = baseShares,
      shareNames = shareNamesMono )

   concResult <- aidsConcav( priceNames = priceNames, totExpName = totExpName,
      data = data, coef = coef, shareNames = shareNames )

   result <- list()
   result$mPercent <- monoResult$mPercent
   result$monotony <- monoResult$monotony
   result$cPercent <- concResult$cPercent
   result$concavity <- concResult$concavity
   result$cMatrices <- concResult$cMatrices

   class( result ) <- "aidsConsist"
   return( result )
}
