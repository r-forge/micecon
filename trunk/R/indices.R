micEconIndex <- function( prices, quantities, base, data, method, na.rm, what ) {

   if( length( prices ) != length( quantities ) ) {
      stop( "arguments 'prices' and 'quantities' must have the same length." )
   }

   checkNames( c( prices, quantities ), names( data ) )

   n <- length( prices )

   numerator <- numeric( nrow( data ) )
   denominator <- numeric( nrow( data ) )

   if( method %in% c( "Laspeyres", "Paasche" ) ) {
      for( i in 1:n ) {
         pt <- data[[ prices[ i ] ]]
         p0 <- mean( data[[ prices[ i ] ]][ base ], na.rm = na.rm )
         qt <- data[[ quantities[ i ] ]]
         q0 <- mean( data[[ quantities[ i ] ]][ base ], na.rm = na.rm )
         if( method == "Laspeyres" ) {
            if( is.na( q0 ) || is.na( p0 ) || all( is.na( pt ) ) ) {
               numerator <- NA
               numerator <- NA
            } else {
               selection <- !is.na( pt )
               numerator[ selection ] <- numerator[ selection ] +
                  pt[ selection ] * q0
               denominator[ selection ] <- denominator[ selection ] + p0 * q0
               if( q0 > 0 ) {
                  numerator[ is.na( pt ) ] <- NA
               }
            }
         } else if( method == "Paasche" ) {
            if( is.na( p0 ) || all( is.na( pt ) ) || all( is.na( qt ) )  ) {
               numerator <- NA
               numerator <- NA
            } else {
               selection <- qt > 0 & !is.na( qt ) & !is.na( pt )
               numerator[ selection ] <- numerator[ selection ] +
                  pt[ selection ] * qt[ selection ]
               denominator[ selection ] <- denominator[ selection ] +
                  p0 * qt[ selection ]
               numerator[ is.na( qt ) ] <- NA
               numerator[ qt > 0  & is.na( pt ) ] <- NA
               denominator[ is.na( qt ) ] <- NA
            }
         }
      }
      result <- numerator / denominator
   } else if( method == "Fisher" ) {
      pL <- priceIndex( prices, quantities, base, data, method = "Laspeyres",
         na.rm = na.rm )
      pP <- priceIndex( prices, quantities, base, data, method = "Paasche",
         na.rm = na.rm )
      result <- sqrt( pL * pP )
   } else {
      stop( paste( "argument 'method' must be either 'Laspeyres', 'Paasche'",
         "or 'Fisher'" ) )
   }
   names( result ) <- rownames( data )
   return( result )
}

priceIndex <- function( prices, quantities, base, data,
   method = "Laspeyres", na.rm = FALSE ) {

   checkNames( c( prices, quantities ), names( data ) )

   result <- micEconIndex( prices, quantities, base, data, method, na.rm,
      "price index" )

   return( result )
}

quantityIndex <- function( prices, quantities, base, data,
   method = "Laspeyres", na.rm = FALSE ) {

   checkNames( c( prices, quantities ), names( data ) )

   result <- micEconIndex( quantities, prices, base, data, method, na.rm,
      "quantity index" )

   return( result )
}


