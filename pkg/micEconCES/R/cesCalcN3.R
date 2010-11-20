cesCalcN3 <- function( xNames, data, coef ) {

   if( coef[ "rho" ] == 0 ) {
      result <- coef[ "gamma_2" ] *
         exp( coef[ "nu" ] * ( coef[ "delta_2" ] *
            ( log( coef[ "gamma_1" ] ) -
               log( coef[ "delta_1" ] * data[[ xNames[ 1 ] ]]^( -coef[ "rho_1" ] ) +
                  ( 1 - coef[ "delta_1" ] ) * data[[ xNames[ 2 ] ]]^( -coef[ "rho_1" ] ) ) /
               coef[ "rho_1" ] ) +
            ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
   } else if( coef[ "rho_1" ] == 0 ) {
      result <- coef[ "gamma_2" ] * 
         ( coef[ "delta_2" ] * coef[ "gamma_1" ]^( -coef[ "rho" ] ) *
            exp( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] )
               )^( -coef[ "rho" ] ) +
            ( 1 - coef[ "delta_2" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] )
   } else {
      result <-
         coef[ "gamma_2" ] * (
            coef[ "delta_2" ] * coef[ "gamma_1" ]^( - coef[ "rho" ] ) *
            ( coef[ "delta_1" ] * data[[ xNames[ 1 ] ]]^( -coef[ "rho_1" ] ) +
               ( 1 - coef[ "delta_1" ] ) * data[[ xNames[ 2 ] ]]^( -coef[ "rho_1" ] )
            )^( coef[ "rho" ] / coef[ "rho_1" ] ) +
            ( 1 - coef[ "delta_2" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] )
         )^( - coef[ "nu" ] / coef[ "rho" ] )
   }

   return( result )
}
