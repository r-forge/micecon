cesCalcN4 <- function( xNames, tName, data, coef ) {

   gammaStar <- coef[ "gamma" ]
   if( "lambda" %in% names( coef ) ) {
      gammaStar <- gammaStar * exp( coef[ "lambda" ] * data[[ tName ]] )
   }

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         if( coef[ "rho_2" ] == 0 ) {
            result <- gammaStar *
               exp( - coef[ "nu" ] *
                  ( coef[ "delta_3" ] * 
                     ( - coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) -
                        ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) ) +
                  ( 1 - coef[ "delta_3" ] ) * 
                     ( - coef[ "delta_2" ] * log( data[[ xNames[ 3 ] ]] ) -
                        ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 4 ] ]] ) ) ) )
         } else {
            result <- gammaStar *
               exp( - coef[ "nu" ] *
                  ( coef[ "delta_3" ] * 
                     ( - coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) -
                        ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) ) +
                  ( 1 - coef[ "delta_3" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
         }
      } else if( coef[ "rho_2" ] == 0 ) {
         result <- gammaStar *
            exp( - coef[ "nu" ] *
               ( coef[ "delta_3" ] * log( B1 ) / coef[ "rho_1" ] +
                  ( 1 - coef[ "delta_3" ] ) *  
                     ( - coef[ "delta_2" ] * log( data[[ xNames[ 3 ] ]] ) -
                        ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 4 ] ]] ) ) ) )
      } else {
         result <- gammaStar *
            exp( - coef[ "nu" ] *
               ( coef[ "delta_3" ] * log( B1 ) / coef[ "rho_1" ] +
                  ( 1 - coef[ "delta_3" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      if( coef[ "rho_2" ] == 0 ) {
         result <- gammaStar *
            ( coef[ "delta_3" ] * 
               exp( coef[ "rho" ] * 
                  ( - coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) -
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) ) ) +
               ( 1 - coef[ "delta_3" ] ) * 
               exp( coef[ "rho" ] * 
                  ( - coef[ "delta_2" ] * log( data[[ xNames[ 3 ] ]] ) -
                     ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 4 ] ]] ) ) ) 
            )^( - ( coef[ "nu" ] / coef[ "rho" ] ) )
      } else {
         result <- gammaStar *
            ( coef[ "delta_3" ] * 
               exp( coef[ "rho" ] * 
                  ( - coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) -
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) ) ) +
               ( 1 - coef[ "delta_3" ] ) * 
                  B2^( coef[ "rho" ] / coef[ "rho_2" ] ) 
            )^( - ( coef[ "nu" ] / coef[ "rho" ] ) )
      }
   } else if( coef[ "rho_2" ] == 0 ) {
      result <- gammaStar *
         ( coef[ "delta_3" ] * B1^( coef[ "rho" ] / coef[ "rho_1" ] ) +
            ( 1 - coef[ "delta_3" ] ) * 
               exp( coef[ "rho" ] * 
                  ( - coef[ "delta_2" ] * log( data[[ xNames[ 3 ] ]] ) -
                     ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 4 ] ]] ) ) ) 
         )^( - ( coef[ "nu" ] / coef[ "rho" ] ) )
   } else {
      result <- gammaStar * ( coef[ "delta_3" ] *
            ( coef[ "delta_1" ] * data[[ xNames[ 1 ] ]]^( -coef[ "rho_1" ] ) +
               ( 1 - coef[ "delta_1" ] ) * data[[ xNames[ 2 ] ]]^( -coef[ "rho_1" ] ) 
            )^( coef[ "rho" ] / coef[ "rho_1" ] ) +
            ( 1 - coef[ "delta_3" ] ) *
            ( coef[ "delta_2" ] * data[[ xNames[ 3 ] ]]^( -coef[ "rho_2" ] ) +
               ( 1 - coef[ "delta_2" ] ) * data[[ xNames[ 4 ] ]]^( -coef[ "rho_2" ] ) 
            )^( coef[ "rho" ] / coef[ "rho_2" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] )
   }

   return( result )
}
