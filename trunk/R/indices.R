priceIndex <- function( prices, quantities, base, data, method = "Laspeyres",
   na.rm = FALSE ) {

   if( length( prices ) != length( quantities ) ) {
      stop( "arguments 'prices' and 'quantities' must have the same length." )
   }

   n <- length( prices )

   numerator <- numeric( nrow( data ) )
   denominator <- numeric( nrow( data ) )

   if( method %in% c( "Laspeyres", "Paasche" ) ) {
      for( i in 1:n ) {
         pt <- with( data, get( prices[ i ] ) )
         p0 <- mean( with( data, get( prices[ i ] ) )[ base ], na.rm = na.rm )
         qt <- with( data, get( quantities[ i ] ) )
         q0 <- mean( with( data, get( quantities[ i ] ) )[ base ], na.rm = na.rm )
         if( method == "Laspeyres" ) {
            if( q0 > 0 ) {
               numerator <- numerator + pt * q0
               denominator <- denominator + p0 * q0
            }
         } else if( method == "Paasche" ) {
            selection <- qt > 0
            numerator[ selection ] <- numerator[ selection ] +
               pt[ selection ] * qt[ selection ]
            denominator[ selection ] <- denominator[ selection ] +
               p0 * qt[ selection ]
         }
      }
      result <- numerator / denominator
   } else if( method == "Fisher" ) {
      pL <- priceIndex( prices, quantities, base, data, method = "Laspeyres" )
      pP <- priceIndex( prices, quantities, base, data, method = "Paasche" )
      result <- sqrt( pL * pP )
   } else {
      stop( paste( "argument 'method' must be either 'Laspeyres', 'Paasche'",
         "or 'Fisher'" ) )
   }
   names( result ) <- rownames( data )
   return( result )
}

quantityIndex <- function( prices, quantities, base, data,
   method = "Laspeyres", na.rm = FALSE ) {

   result <- priceIndex( quantities, prices, base, data, method, na.rm )

   return( result )
}


