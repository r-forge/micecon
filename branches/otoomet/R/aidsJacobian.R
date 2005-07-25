aidsJacobian <- function( allCoef, pNames, xtName, data = NULL,
      omitLast = TRUE, alpha0 = 0 ) {
   nObs <- nrow( data )
   coef <- aidsCoef( allCoef )
   nGoods <- length( coef$alpha )
   hom <- all.equal( rowSums( coef$gamma ), rep( 0, nGoods ) ) == TRUE
   sym <- all.equal( coef$gamma, t( coef$gamma ) ) == TRUE
   lnp <- aidsPx( "TL", pNames, coef = coef, data = data, alpha0 = alpha0 )
   result <- matrix( 0, nrow = nObs * ( nGoods - 1 ),
      ncol = ( nGoods + 2 ) * ( nGoods - 1 ) )
   for( eq in 1:( nGoods - 1 ) ) {
      myRows <- ( ( eq - 1 ) * nObs + 1 ):( eq * nObs )
      # derivatives of alphas
      for( i in 1:( nGoods - 1 ) ) {
         myCol <- ( i - 1 ) * ( nGoods + 2 ) + 1
         result[ myRows, myCol ] <- ( i == eq ) -
            coef$beta[ eq ] *
            ( log( data[[ pNames[ i ] ]] ) -
            log( data[[ pNames[ nGoods ] ]] ) )
      }
      # derivatives of betas
      myCol <- ( eq - 1 ) * ( nGoods + 2 ) + 2
      result[ myRows, myCol ] <- log( data[[ xtName ]] ) - lnp
      # derivatives of gammas
      for( i in 1:( nGoods - 1 ) ) {
         for( j in 1:( nGoods - hom ) ) {
            myCol <- ( i - 1 ) * ( nGoods + 2 ) + 2 + j
            result[ myRows, myCol ] <-
               ( i == eq ) * ( log( data[[ pNames[ j ] ]] ) -
                  hom * log( data[[ pNames[ nGoods ] ]] ) ) -
               0.5 * coef$beta[ eq ] *
               ( log( data[[ pNames[ i ] ]] ) -
                  log( data[[ pNames[ nGoods ] ]] ) ) *
               ( log( data[[ pNames[ j ] ]] ) -
                  hom * log( data[[ pNames[ nGoods ] ]] ) )
         }
      }
   }
   delCols <- NULL
   for( i in 1:( nGoods - 1 ) ) {
      if( hom ) {
         delCols <- c( delCols, i * ( nGoods + 2 ) )
      }
      if( sym && i >= 2 ) {
         for( j in 1:( i - 1 ) ) {
            delCol <- ( i - 1 ) * ( nGoods + 2 ) + 2 + j
            stayCol <- ( j - 1 ) * ( nGoods + 2 ) + 2 + i
            result[ , stayCol ] <- result[ , stayCol ] + result[ , delCol ]
            delCols <- c( delCols, delCol )
         }
      }
   }
   if( !is.null( delCols ) ) {
      result <- result[ , -delCols ]
   }
   return( result )
}


if( FALSE ) {

function( nGoods, data, b, alpha0 = 0, omitLast = TRUE ) {
   on.exit( print( result ) )
   nObs <- nrow( data )
   coef <- aidsCoef( b )
   for( i in 1:nGoods ) {
      data[[ paste( "p", i, sep = "" ) ]] <-
         exp( data[[ paste( "lp", i, sep = "" ) ]] )
   }
   data$lnp <- aidsPx( "TL", paste( "p", c( 1:nGoods ), sep = "" ),
      data = data, alpha0 = alpha0, coef = coef )
   result <- matrix( 0, nrow = nObs * ( nGoods - omitLast ),
      ncol = ( nGoods - omitLast ) * ( nGoods + 2 ) )
   for( i in 1:( nGoods - omitLast ) ) {
      myRows <- ( ( i - 1 ) * nObs + 1 ):( i * nObs )
      # derivatives of alphas
      result[ myRows, ( i - 1 ) * ( nGoods + 2 ) + 1 ] <- 1
      for( j in 1:( nGoods - omitLast ) ) {
         result[ myRows, ( j - 1 ) * ( nGoods + 2 ) + 1 ] <-
            result[ myRows, ( j - 1 ) * ( nGoods + 2 ) + 1 ] -
            coef$beta[ i ] * data[[ paste( "lp", j, sep = "" ) ]]
      }
      # derivatives of betas
      result[ myRows, ( i - 1 ) * ( nGoods + 2 ) + 2 ] <-
         log( data$xt ) - data$lnp
      # derivatives of gammas
      for( k in 1:( nGoods - omitLast ) ) {
         result[ myRows, ( i - 1 ) * ( nGoods + 2 ) + 2 + k ] <-
            data[[ paste( "lp", k, sep = "" ) ]]
         for( j in 1:( nGoods - omitLast ) ) {
            result[ myRows, ( j - 1 ) * ( nGoods + 2 ) + 2 + k ] <-
               result[ myRows, ( j - 1 ) * ( nGoods + 2 ) + 2 + k ] -
               0.5 * coef$beta[ i ] * data[[ paste( "lp", j, sep = "" ) ]] *
               data[[ paste( "lp", k, sep = "" ) ]]
         }
      }
   }
   on.exit( )
   return( result )
}
}

if( FALSE ) {

   sysData <- data.frame( xt = data[[ xtName ]],
      lxtr = ( log( data[[ xtName ]] ) - lnp ) )


      # calculating log of "real" (deflated) total expenditure
      sysData$lxtr <- log( data[[ xtName ]] ) -
         aidsPx( "TL", pNames, wNames, data = data,
         alpha0 = alpha0, coef = aidsCoef( est$b ) )
      # calculating matrix G
      Gmat <- cbind( rep( 1, nObs ), sysData$lxtr )
      for( i in 1:( nGoods ) ) {
         Gmat <- cbind( Gmat, sysData[[ paste( "lp", i, sep = "" ) ]] )
      }
      # testing matrix G
      if( FALSE ) {
         for( i in 1:( nGoods - 1 ) ) {
            print( est$eq[[ i ]]$fitted - Gmat %*%
               est$b[ ( ( i - 1 ) * ( nGoods + 2 ) + 1 ):( i * (nGoods + 2 ) ) ] )
         }
      }
      # calculating matrix J
      jacobian <- aidsJacobian( nGoods, sysData, est$b, alpha0 = alpha0 )
      Jmat <- ( diag( nGoods - 1 ) %x% t( Gmat ) ) %*% jacobian
      JmatInv <- solve( Jmat )
      bcov <- JmatInv %*% ( est$rcov %x% ( t( Gmat ) %*% Gmat ) ) %*%
         t( JmatInv )
}