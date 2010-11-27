cesCalcN4 <- function( xNames, data, coef ) {

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      result <- coef[ "gamma" ] *
         exp( - coef[ "nu" ] *
            ( coef[ "delta_3" ] * log( B1 ) / coef[ "rho_1" ] +
               ( 1 - coef[ "delta_3" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
   } else {
      result <- coef[ "gamma" ] * ( coef[ "delta_3" ] *
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
