aidsElas <- function( coef, shares, prices = NULL, formula = "AIDS",
   quantNames = NULL, priceNames = NULL, coefVcov = NULL, df = NULL ) {

   nGoods <- length( coef$alpha )

   if( length( coef$alpha ) != length( coef$beta ) ) {
      stop( "arguments 'alpha' and 'beta' must have the same length" )
   } else if( nrow( coef$gamma ) != ncol( coef$gamma ) ) {
      stop( "argument 'gamma' must be a square matrix" )
   } else if( length( coef$alpha ) != nrow( coef$gamma ) ) {
      stop( "number of rows of argument 'gamma' must be equal",
         " to the length of argument 'alpha'" )
   } else if(  length( coef$alpha ) != length( shares ) ) {
      stop( "arguments 'alpha' and 'shares' must have the same length" )
   } else if(  length( coef$alpha ) != length( prices ) && !is.null( prices ) ) {
      stop( "arguments 'alpha' and 'prices' must have the same length" )
   }
   if( is.null( quantNames ) ) {
      quantNames <- .aidsQuantNames( shares, coef, nGoods )
   } else {
      if( length( quantNames ) != nGoods ) {
         stop( "argument 'quantNames' must have ", nGoods, " elements" )
      }
   }
   if( is.null( priceNames ) ) {
      priceNames <- .aidsPriceNames( prices, coef, nGoods )
   } else {
      if( length( priceNames ) != nGoods ) {
         stop( "argument 'priceNames' must have ", nGoods, " elements" )
      }
   }


   if( formula %in% c( "AIDS" ) ) {
      if( is.null( prices ) ) {
         stop( "the 'AIDS' formula requires argument 'prices'" )
      }
   } else if( formula %in% c( "Ch", "EU" ) ) {
      if( !is.null( prices ) ) {
         warning( "the 'Ch' and 'EU' formulas do not require argument 'prices'" )
      }
   }

   ela <- list()
   ela$formula <- formula

   ones <- rep( 1, nGoods )

   if( formula == "AIDS" ) {
      ela$exp <- ones + coef$beta/shares
      ela$hicks <- -diag( 1, nGoods, nGoods ) +
         ones %*% t( shares ) +
         coef$gamma / ( shares %*% t( ones ) ) -
         coef$beta %*% t( ones ) *
         ( ones %*% t( coef$alpha ) -
         ones %*% t( shares )+
         ones %*% t( coef$gamma %*% log( prices ))) /
         ( shares %*% t( ones ) )
      ela$marshall <- ela$hicks - ( ela$exp %*% t( ones ) ) *
         ( ones %*% t( shares ))
   } else if(formula=="Ch") {
      ela$exp <- ones + coef$beta / shares
      ela$hicks <- -diag( 1, nGoods, nGoods ) + coef$gamma /
         ( shares %*% t( ones ) ) +
         ones %*% t( shares )
      ela$marshall <- ela$hicks - ( ela$exp %*% t( ones ) ) *
         ( ones %*% t(shares))
   } else if( formula == "EU" ) {
      ela$exp <- ones + coef$beta / shares
      ela$marshall <- -diag( 1, nGoods, nGoods ) + coef$gamma /
         ( shares %*% t( ones ) )
      ela$hicks <- ela$marshall + ( ela$exp %*% t( ones ) ) *
         ( ones %*% t(shares))
   } else if( formula == "B1 not implemented" ) {
      ela$exp <- array( 1, c( nGoods ) ) + coef$beta / shares
      ela$marshall <- matrix( NA, nGoods, nGoods )
      for( i in 1:nGoods ) {
         for( j in 1:nGoods ) {
            ela$marshall[ i, j ] <- -( i == j ) + coef$gamma[ i, j ] / shares[ i ] -
               coef$beta[ i ] * shares[ j ] / shares[ i ]
            for( k in 1:nGoods ) {
               ela$marshall[ i, j ] <- ela$marshall[ i, j ] - coef$beta[ i ] *
                  shares[ k ] * log( prices[ k ] ) * ( 0 + ( k == j ) ) / shares[ i ]
            }
         }
      }
      ela$hicks <- ela$marshall + ( ela$exp %*% t( array( 1, c( nGoods ) ) ) ) *
         ( array( 1, c( nGoods )) %*% t(shares))
   } else {
      stop( "formula '", as.character( formula ), "' is not supported" )
   }
   names( ela$exp )         <- quantNames
   rownames( ela$hicks )    <- quantNames
   colnames( ela$hicks )    <- priceNames
   rownames( ela$marshall ) <- quantNames
   colnames( ela$marshall ) <- priceNames
   if( !is.null( coefVcov ) && formula %in% c( "AIDS" ) ) {
      jacobian <- aidsElasJacobian( coef = coef, share = shares, price = prices,
         formula = formula, quantNames = quantNames, priceNames = priceNames )
      ela$allVcov      <- jacobian$all      %*% coefVcov %*% t( jacobian$all )
      ela$expVcov      <- jacobian$exp      %*% coefVcov %*% t( jacobian$exp )
      ela$hicksVcov    <- jacobian$hicks    %*% coefVcov %*% t( jacobian$hicks )
      ela$marshallVcov <- jacobian$marshall %*% coefVcov %*% t( jacobian$marshall )
      # standard errors
      ela$expStEr      <- diag( ela$expVcov )^0.5
      ela$hicksStEr    <- matrix( diag( ela$hicksVcov )^0.5,
         ncol = nGoods, byrow = TRUE )
      ela$marshallStEr <-  matrix( diag( ela$marshallVcov )^0.5,
         ncol = nGoods, byrow = TRUE )
      # dim names for standard errors
      names( ela$expStEr )         <- names( ela$exp )
      dimnames( ela$hicksStEr )    <- dimnames( ela$hicks )
      dimnames( ela$marshallStEr ) <- dimnames( ela$marshall )
      # t-values
      ela$expTval      <- ela$exp      / ela$expStEr
      ela$hicksTval    <- ela$hicks    / ela$hicksStEr
      ela$marshallTval <- ela$marshall / ela$marshallStEr
      if( !is.null( df ) ) {
         ela$expPval <- 2 * pt( abs( ela$expTval ), df,
            lower.tail = FALSE )
         ela$hicksPval <- 2 * pt( abs( ela$hicksTval ), df,
            lower.tail = FALSE )
         ela$marshallPval <- 2 * pt( abs( ela$marshallTval ), df,
            lower.tail = FALSE )
      }
   }
   class( ela ) <- "aidsElas"
   return( ela )
}
