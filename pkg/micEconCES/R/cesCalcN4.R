cesCalcN4 <- function( xNames, data, coef ) {

      result <- coef[ "gamma" ] * ( coef[ "delta_3" ] *
            ( coef[ "delta_1" ] * data[[ xNames[ 1 ] ]]^( -coef[ "rho_1" ] ) +
               ( 1 - coef[ "delta_1" ] ) * data[[ xNames[ 2 ] ]]^( -coef[ "rho_1" ] ) 
            )^( coef[ "rho" ] / coef[ "rho_1" ] ) +
            ( 1 - coef[ "delta_3" ] ) *
            ( coef[ "delta_2" ] * data[[ xNames[ 3 ] ]]^( -coef[ "rho_2" ] ) +
               ( 1 - coef[ "delta_2" ] ) * data[[ xNames[ 4 ] ]]^( -coef[ "rho_2" ] ) 
            )^( coef[ "rho" ] / coef[ "rho_2" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] )

   return( result )
}
