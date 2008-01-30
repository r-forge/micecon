.aidsCheckCoef <- function( coef, argCoef = "coef",
      nGoods = NA, argGoods = NULL,
      nShifters = NA, argShifters = NULL ) {

   # checking arguments of nGoods and argGoods of *this* function
   if( any( !is.na( nGoods ) ) ){
      if( is.null( argGoods ) ){
         stop( "internal error: 'nGoods' is specified",
            " but 'argGoods' is not available" )
      }
      if( length( nGoods ) != length( argGoods ) ){
         stop( "internal error: 'nGoods' and 'argGoods'",
            " have not the same length" )
      }
      for( i in 1:length( nGoods ) ){
         if( !is.numeric( nGoods[ i ] ) || length( nGoods[ i ] ) != 1 ){
            stop( "internal error: 'nGoods['", i , "] must be",
               " a numeric scalar" )
         }
         if( !is.character( argGoods[ i ] ) ){
            stop( "internal error: 'argGoods['", i , "] must be",
               " a character string" )
         }
      }
   }

  # checking arguments of nShifters and argShifters of *this* function
   if( any( !is.na( nShifters ) ) ){
      if( is.null( argShifters ) ){
         stop( "internal error: 'nShifters' is specified",
            " but 'argShifters' is not available" )
      }
      if( length( nShifters ) != length( argShifters ) ){
         stop( "internal error: 'nShifters' and 'argShifters'",
            " have not the same length" )
      }
      for( i in 1:length( nShifters ) ){
         if( !is.numeric( nShifters[ i ] ) || length( nShifters[ i ] ) != 1 ){
            stop( "internal error: 'nShifters['", i , "] must be",
               " a numeric scalar" )
         }
         if( !is.character( argShifters[ i ] ) ){
            stop( "internal error: 'argShifters['", i , "] must be",
               " a character string" )
         }
      }
   }

   ## checking coefficients
   # alpha
   if( is.null( coef$alpha ) ){
      return( paste( "'", argCoef, "$alpha' is missing", sep = "" ) )
   }
   if( !is.numeric( coef$alpha ) ){
      return( paste( "'", argCoef, "$alpha' must be numeric vector", sep = "" ) )
   }

   # beta
   if( is.null( coef$beta ) ){
      return( paste( "'", argCoef, "$beta' is missing", sep = "" ) )
   }
   if( !is.numeric( coef$beta ) ){
      return( paste( "'", argCoef, "$beta' must be numeric vector", sep = "" ) )
   }
   if( length( coef$alpha ) != length( coef$beta ) ) {
      return( paste( "'", argCoef, "$alpha' and '", argCoef, "$beta'",
         " must have the same length", sep = "" ) )
   }

   # gamma
   if( is.null( coef$gamma ) ){
      return( paste( "'", argCoef, "$gamma' is missing", sep = "" ) )
   }
   if( !is.matrix( coef$gamma ) ){
      return( paste( "'", argCoef, "$gamma' must be a matrix", sep = "" ) )
   }
   if( nrow( coef$gamma ) != ncol( coef$gamma ) ) {
      return( paste( "argument '", argCoef, "$gamma' must be a square matrix",
         sep = "" ) )
   }
   if( length( coef$alpha ) != nrow( coef$gamma ) ) {
      return( paste( "the number of rows of '", argCoef, "$gamma'",
         " must be equal to the length of '", argCoef, "$alpha'",
         sep = "" ) )
   }

   # delta
   if( !is.null( coef$delta ) ){
      if( !is.matrix( coef$delta ) ){
         return( paste( "'", argCoef, "$delta' must be a matrix", sep = "" ) )
      }
      if( length( coef$alpha ) != nrow( coef$delta ) ) {
         return( paste( "the number of rows of '", argCoef, "$delta'",
            " must be equal to the length of '", argCoef, "$alpha'",
            sep = "" ) )
      }
   }

   # checking Goods
   if( any( !is.na( nGoods ) ) ){
      for( i in 1:length( nGoods ) ){
         if( nGoods[i] != length( coef$alpha ) && !is.na( nGoods[i] ) ) {
            return( paste( "'", argCoef, "$alpha' and '", argGoods[i],
               "' must have the same length", sep = "" ) )
         }
      }
   }

   # checking demand shifters
   if( any( !is.na( nShifters ) ) ){
      if( is.null( coef$delta ) ){
         stop( "'", argCoef, "$delta' is not available" )
      }
      for( i in 1:length( nGoods ) ){
         if( nShifters != ncol( coef$delta ) && !is.na( nShifters[i] ) ) {
            return( paste( "the number of columns of '", argCoef, "$delta'",
               " must be equal to the length of ", argShifters, sep = "" ) )
         }
      }
   }

   return( NULL )
}