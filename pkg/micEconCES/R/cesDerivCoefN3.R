# derivatives with respect to gamma_1
cesDerivCoefN3Gamma1 <- function( coef, B1, B ) {

   result <- coef[ "gamma_2" ] * coef[ "nu" ] * 
      coef[ "delta_2" ] * coef[ "gamma_1" ]^( - coef[ "rho" ] - 1 ) * 
      B1^( coef[ "rho" ] / coef[ "rho_1" ] ) * 
      B^( - coef[ "nu" ] / coef[ "rho" ] - 1 )

   return( result )
}


# derivatives with respect to gamma_2
cesDerivCoefN3Gamma2 <- function( coef, B ) {

   result <- B^( - coef[ "nu" ] / coef[ "rho" ] )

   return( result )
}


# derivatives with respect to delta_1
cesDerivCoefN3Delta1 <- function( coef, B1, B, data, xNames ) {

   result <- - coef[ "gamma_2" ] * coef[ "nu" ] * coef[ "delta_2" ] * 
      coef[ "gamma_1" ]^( -coef[ "rho" ] ) * 
      B1^( coef[ "rho" ] / coef[ "rho_1" ] - 1 ) / coef[ "rho_1" ] * 
      B^( -coef[ "nu" ] / coef[ "rho" ] - 1 ) *
      ( data[[ xNames[ 1 ] ]]^( -coef[ "rho_1" ] ) 
         - data[[ xNames[ 2 ] ]]^( -coef[ "rho_1" ] ) )

   return( result )
}


# derivatives with respect to delta_2
cesDerivCoefN3Delta2 <- function( coef, B1, B, data, xNames ) {

   result <- -coef[ "gamma_2" ] * ( coef[ "nu" ] / coef[ "rho" ] ) *
      ( coef[ "gamma_1" ]^( -coef[ "rho" ] ) * 
      B1^( coef[ "rho" ] / coef[ "rho_1" ] ) - 
      data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) ) *
      B^( -coef[ "nu" ] / coef[ "rho" ] - 1 )

   return( result )
}


# derivatives with respect to rho_1
cesDerivCoefN3Rho1 <- function( coef, B1, B, data, xNames ) {

   result <- -coef[ "gamma_2" ] * coef[ "nu" ] * coef[ "delta_2" ] / 
      coef[ "rho_1" ] * coef[ "gamma_1" ]^( -coef[ "rho" ] ) * 
      B1^( coef[ "rho" ] / coef[ "rho_1" ] ) * 
      B^( -coef[ "nu" ] / coef[ "rho" ] - 1 ) *
      ( - log( B1 ) / coef[ "rho_1" ] +
      ( -coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) * 
         data[[ xNames[ 1 ] ]]^( -coef[ "rho_1" ] ) -
         ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) *
         data[[ xNames[ 2 ] ]]^( -coef[ "rho_1" ] ) ) / B1 )

   return( result )
}


# derivatives with respect to rho
cesDerivCoefN3Rho <- function( coef, B1, B, data, xNames ) {

   result <- -coef[ "gamma_2" ] * ( coef[ "nu" ] / coef[ "rho" ] ) * 
      B^( -coef[ "nu" ] / coef[ "rho" ] ) *
      ( -log( B ) / coef[ "rho" ] +
         B^(-1) * ( coef[ "delta_2" ] * coef[ "gamma_1" ]^( -coef[ "rho" ] ) * 
         B1^( coef[ "rho" ] / coef[ "rho_1" ] ) *
            ( - log( coef[ "gamma_1" ] ) + log( B1 ) / coef[ "rho_1" ] ) -
         ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 3 ] ]] ) *
         data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) ) )

   return( result )
}


# derivatives with respect to nu
cesDerivCoefN3Nu <- function( coef, B ) {

   result <- - coef[ "gamma_2" ] * log( B ) * 
      B^( -coef[ "nu" ] / coef[ "rho" ] ) / coef[ "rho" ]

   return( result )
}

