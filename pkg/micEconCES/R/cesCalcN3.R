cesCalcN3 <- function( xNames, data, coef ) {

   result <-
      coef[ "gamma_2" ] * (
         coef[ "delta_2" ] * coef[ "gamma_1" ]^( - coef[ "rho" ] ) *
         ( coef[ "delta_1" ] * data[[ xNames[ 1 ] ]]^( -coef[ "rho_1" ] ) +
            ( 1 - coef[ "delta_1" ] ) * data[[ xNames[ 2 ] ]]^( -coef[ "rho_1" ] )
         )^( coef[ "rho" ] / coef[ "rho_1" ] ) +
         ( 1 - coef[ "delta_2" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] )
      )^( - coef[ "nu" ] / coef[ "rho" ] )

   return( result )
}
