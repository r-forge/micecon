aidsConsist <- function( priceNames, totExpName, coef, data,
      priceIndex = "TL", basePrices = NULL, baseShares = NULL,
      shareNames = NULL ) {

   result <- list()

   result$addingUp <- all(
      all.equal( sum( coef$alpha ), 1 ) == TRUE,
      all.equal( sum( coef$beta ), 0 ) == TRUE,
      all.equal( colSums( coef$gamma ), rep( 0, 4 ),
         check.attributes = FALSE ) == TRUE
      )

   result$homogeneity <- all.equal( rowSums( coef$gamma ),
      rep( 0, 4 ), check.attributes = FALSE ) == TRUE

   result$symmetry <- isSymmetric( coef$gamma, tol = 1e-10,
      check.attributes = FALSE ) == TRUE

   result$mono <- aidsMono( priceNames = priceNames, totExpName = totExpName,
      data = data, coef = coef, priceIndex = priceIndex,
      basePrices = basePrices, baseShares = baseShares )

   if( priceIndex == "TL" ) {
      result$concav <- aidsConcav( priceNames = priceNames,
         totExpName = totExpName, data = data, coef = coef,
         shareNames = shareNames )
   }

   class( result ) <- "aidsConsist"
   return( result )
}
