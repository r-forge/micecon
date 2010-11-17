# calculate part "B1"
cesDerivCoefN3B1 <- function( coef, data, xNames ) {

   result <- coef[ "delta_1" ] * data[[ xNames[ 1 ] ]]^( -coef[ "rho_1" ] ) +
      ( 1 - coef[ "delta_1" ] ) * data[[ xNames[ 2 ] ]]^( -coef[ "rho_1" ] )

   return( result )
}


# calculate part "B"
cesDerivCoefN3B <- function( coef, data, xNames ) {

   B1 <- cesDerivCoefN3B1( coef = coef, data = data, xNames = xNames )

   result <- coef[ "delta_2" ] * coef[ "gamma_1" ]^( -coef[ "rho" ] ) * 
      B1^( coef[ "rho" ] / coef[ "rho_1" ] ) + 
      ( 1 - coef[ "delta_2" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] )

   return( result )
}

# derivatives with respect to gamma_1
cesDerivCoefN3Gamma1 <- function( coef, data, xNames ) {

   B1 <- cesDerivCoefN3B1( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN3B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         result <- ( coef[ "gamma_2" ] / coef[ "gamma_1" ] ) *
            coef[ "nu" ] * coef[ "delta_2" ] *
            exp( coef[ "nu" ] * ( coef[ "delta_2" ] *
                  ( log( coef[ "gamma_1" ] ) + 
                  coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) ) +
               ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      } else {
         result <- ( coef[ "gamma_2" ] / coef[ "gamma_1" ] ) * 
            coef[ "nu" ] * coef[ "delta_2" ] *
               exp( - coef[ "nu" ] * 
                  ( - coef[ "delta_2" ] * 
                     ( log( coef[ "gamma_1" ] ) - log( B1 ) / coef[ "rho_1" ] ) -
                     ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      result <- coef[ "gamma_2" ] * coef[ "nu" ] * coef[ "delta_2" ] * 
         coef[ "gamma_1" ]^( -coef[ "rho" ] - 1 ) *
         exp( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
            ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) 
         )^( -coef[ "rho" ] ) *
         ( coef[ "delta_2" ] * coef[ "gamma_1" ]^( -coef[ "rho" ] ) *
            exp( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
               ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) 
            )^( -coef[ "rho" ] ) +
            ( 1 - coef[ "delta_2" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] - 1 )
   } else {
      result <- coef[ "gamma_2" ] * coef[ "nu" ] * 
         coef[ "delta_2" ] * coef[ "gamma_1" ]^( - coef[ "rho" ] - 1 ) * 
         B1^( coef[ "rho" ] / coef[ "rho_1" ] ) * 
         B^( - coef[ "nu" ] / coef[ "rho" ] - 1 )
   }

   return( result )
}


# derivatives with respect to gamma_2
cesDerivCoefN3Gamma2 <- function( coef, data, xNames ) {

   B1 <- cesDerivCoefN3B1( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN3B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         result <- exp( coef[ "nu" ] * 
            ( coef[ "delta_2" ] *
               ( log( coef[ "gamma_1" ] ) + 
                  coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) ) +
               ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      } else {
         result <-
            exp( coef[ "nu" ] * 
               ( coef[ "delta_2" ] * 
                  ( log( coef[ "gamma_1" ] ) - log( B1 ) / coef[ "rho_1" ] ) +
                  ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      result <-
         ( coef[ "delta_2" ] * coef[ "gamma_1" ]^( -coef[ "rho" ] ) *
            exp( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
               ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) )^( -coef[ "rho" ] ) +
            ( 1 - coef[ "delta_2" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] )
   } else {
      result <- B^( - coef[ "nu" ] / coef[ "rho" ] )
   }
   return( result )
}


# derivatives with respect to delta_1
cesDerivCoefN3Delta1 <- function( coef, data, xNames ) {

   B1 <- cesDerivCoefN3B1( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN3B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         result <- - coef[ "gamma_2" ] * coef[ "nu" ] * coef[ "delta_2" ] *
               ( - log( data[[ xNames[ 1 ] ]] ) + log( data[[ xNames[ 2 ] ]] ) ) *
               exp( - coef[ "nu" ] * 
                  ( - coef[ "delta_2" ] *
                     ( log( coef[ "gamma_1" ] ) + 
                        coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
                        ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) 
                     ) -
                     ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 3 ] ]] ) 
                  ) )
      } else {
         result <- - coef[ "gamma_2" ] * coef[ "nu" ] * coef[ "delta_2" ] *
            ( data[[ xNames[ 1 ] ]]^(-coef[ "rho_1" ]) - 
               data[[ xNames[ 2 ] ]]^(-coef[ "rho_1" ]) ) /
            ( coef[ "rho_1" ] * B1 ) *
            exp( - coef[ "nu" ] * ( - coef[ "delta_2" ] * 
                  ( log( coef[ "gamma_1" ] ) - log( B1 ) / coef[ "rho_1" ] ) -
               ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      result <-
         coef[ "gamma_2" ] * coef[ "nu" ] * coef[ "delta_2" ] * 
         coef[ "gamma_1" ]^( -coef[ "rho" ] ) *
         exp( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
            ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) 
         )^( -coef[ "rho" ] ) *
         ( log( data[[ xNames[ 1 ] ]] ) - log( data[[ xNames[ 2 ] ]] ) ) *
         ( coef[ "delta_2" ] * coef[ "gamma_1" ]^( -coef[ "rho" ] ) *
            exp( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
               ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) 
            )^( -coef[ "rho" ] ) +
            ( 1 - coef[ "delta_2" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] - 1 )
   } else {
      result <- - coef[ "gamma_2" ] * coef[ "nu" ] * coef[ "delta_2" ] * 
         coef[ "gamma_1" ]^( -coef[ "rho" ] ) * 
         B1^( coef[ "rho" ] / coef[ "rho_1" ] - 1 ) / coef[ "rho_1" ] * 
         B^( -coef[ "nu" ] / coef[ "rho" ] - 1 ) *
         ( data[[ xNames[ 1 ] ]]^( -coef[ "rho_1" ] ) 
            - data[[ xNames[ 2 ] ]]^( -coef[ "rho_1" ] ) )
   }

   return( result )
}


# derivatives with respect to delta_2
cesDerivCoefN3Delta2 <- function( coef, data, xNames ) {

   B1 <- cesDerivCoefN3B1( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN3B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         result <- - coef[ "gamma_2" ] * coef[ "nu" ] * 
            ( - log( coef[ "gamma_1" ] ) -
               coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) -
               ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) +
               log( data[[ xNames[ 3 ] ]] ) ) *
            exp( coef[ "nu" ] * ( coef[ "delta_2" ] *
               ( log( coef[ "gamma_1" ] ) + 
                  coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) ) +
               ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      } else {
         result <- - coef[ "gamma_2" ] * coef[ "nu" ] *
            ( - log( coef[ "gamma_1" ] ) + log( B1 ) / coef[ "rho_1" ] + 
               log( data[[ xNames[ 3 ] ]] ) ) *
            exp( coef[ "nu" ] * ( coef[ "delta_2" ] * 
                  ( log( coef[ "gamma_1" ] ) - log( B1 ) / coef[ "rho_1" ] ) +
               ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      result <-
         - coef[ "gamma_2" ] * ( coef[ "nu" ] / coef[ "rho" ] ) *
         ( coef[ "gamma_1" ]^( -coef[ "rho" ] ) * 
            exp( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
               ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) 
            )^( -coef[ "rho" ] ) -
            data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) ) *
         ( coef[ "delta_2" ] * coef[ "gamma_1" ]^( -coef[ "rho" ] ) *
            exp( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
               ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) 
            )^( -coef[ "rho" ] ) +
            ( 1 - coef[ "delta_2" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] - 1 )
   } else {
      result <- -coef[ "gamma_2" ] * ( coef[ "nu" ] / coef[ "rho" ] ) *
         ( coef[ "gamma_1" ]^( -coef[ "rho" ] ) * 
         B1^( coef[ "rho" ] / coef[ "rho_1" ] ) - 
         data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) ) *
         B^( -coef[ "nu" ] / coef[ "rho" ] - 1 )
   }

   return( result )
}


# derivatives with respect to rho_1
cesDerivCoefN3Rho1 <- function( coef, data, xNames ) {

   B1 <- cesDerivCoefN3B1( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN3B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         result <- - coef[ "gamma_2" ] * coef[ "nu" ] * coef[ "delta_2" ] * 
            0.5 * ( coef[ "delta_1" ] * ( log( data[[ xNames[ 1 ] ]] ) )^2 +
               ( 1 - coef[ "delta_1" ] ) * ( log( data[[ xNames[ 2 ] ]] ) )^2 -
               ( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) )^2 ) *
            exp( coef[ "nu" ] * ( coef[ "delta_2" ] *
               ( log( coef[ "gamma_1" ] ) +
                  coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) ) +
               ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      } else {
         result <- - coef[ "gamma_2" ] * coef[ "nu" ] * 
            ( coef[ "delta_2" ] / coef[ "rho_1" ] ) *
            ( - log( B1 ) / coef[ "rho_1" ] + 
               ( - coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) * 
                  data[[ xNames[ 1 ] ]]^(-coef[ "rho_1" ]) -
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) *
                  data[[ xNames[ 2 ] ]]^(-coef[ "rho_1" ]) ) / B1 ) *
            exp( - coef[ "nu" ] * ( - coef[ "delta_2" ] * 
                  ( log( coef[ "gamma_1" ] ) - log( B1 ) / coef[ "rho_1" ] ) -
               ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      result <- - coef[ "gamma_2" ] * coef[ "nu" ] * coef[ "delta_2" ] * 
         coef[ "gamma_1" ]^(-coef[ "rho" ]) *
         exp( -coef[ "rho" ] * 
            ( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
               ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) ) ) *
         0.5 * ( ( coef[ "delta_1" ] * ( log( data[[ xNames[ 1 ] ]] ) )^2 +
               ( 1 - coef[ "delta_1" ] ) * ( log( data[[ xNames[ 2 ] ]] ) )^2 ) -
            ( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
               ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) )^2 ) *
         ( coef[ "delta_2" ] * coef[ "gamma_1" ]^(-coef[ "rho" ]) *
            exp( -coef[ "rho" ] * 
               ( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) ) ) +
            ( 1 - coef[ "delta_2" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) 
         )^( -coef[ "nu" ]/coef[ "rho" ] - 1 )
   } else {
      result <- -coef[ "gamma_2" ] * coef[ "nu" ] * coef[ "delta_2" ] / 
         coef[ "rho_1" ] * coef[ "gamma_1" ]^( -coef[ "rho" ] ) * 
         B1^( coef[ "rho" ] / coef[ "rho_1" ] ) * 
         B^( -coef[ "nu" ] / coef[ "rho" ] - 1 ) *
         ( - log( B1 ) / coef[ "rho_1" ] +
         ( -coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) * 
            data[[ xNames[ 1 ] ]]^( -coef[ "rho_1" ] ) -
            ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) *
            data[[ xNames[ 2 ] ]]^( -coef[ "rho_1" ] ) ) / B1 )
   }

   return( result )
}


# derivatives with respect to rho
cesDerivCoefN3Rho <- function( coef, data, xNames ) {

   B1 <- cesDerivCoefN3B1( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN3B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         result <- coef[ "gamma_2" ] * coef[ "nu" ] *
            ( -0.5 * 
               ( coef[ "delta_2" ] *
                  ( - log( coef[ "gamma_1" ] ) -
                     coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) -
                     ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) 
                  )^2 +
                  ( 1- coef[ "delta_2" ] ) * ( log( data[[ xNames[ 3 ] ]] ) )^2 
               ) +
               0.5 * 
               ( coef[ "delta_2" ] *
                  ( - log( coef[ "gamma_1" ] ) -
                     coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) -
                     ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) ) -
                  ( 1- coef[ "delta_2" ] ) * log( data[[ xNames[ 3 ] ]] ) )^2 
            ) *
            exp( coef[ "nu" ] * ( - coef[ "delta_2" ] *
               ( - log( coef[ "gamma_1" ] ) - 
                  coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) -
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) ) +
               ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      } else {
         result <- - coef[ "gamma_2" ] * coef[ "nu" ] *
            ( 0.5 * ( coef[ "delta_2" ] * 
                  ( - log( coef[ "gamma_1" ] ) + log( B1 ) / coef[ "rho_1" ] )^2 +
                  ( 1 - coef[ "delta_2" ] ) * ( log( data[[ xNames[ 3 ] ]] ) )^2 ) -
               0.5 * ( coef[ "delta_2" ] * 
                  ( - log( coef[ "gamma_1" ] ) + log( B1 ) / coef[ "rho_1" ] ) -
               ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 3 ] ]] ) )^2 ) *
            exp( coef[ "nu" ] * ( - coef[ "delta_2" ] * 
                  ( - log( coef[ "gamma_1" ] ) + log( B1 ) / coef[ "rho_1" ] ) +
               ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      result <- coef[ "gamma_2" ] * ( 
         ( coef[ "nu" ] / coef[ "rho" ]^2 ) * 
         log( coef[ "delta_2" ] * coef[ "gamma_1" ]^( -coef[ "rho" ] ) *
            ( exp( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) ) 
            )^( -coef[ "rho" ] ) +
            ( 1 - coef[ "delta_2" ] ) * 
            data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) ) *
         ( coef[ "delta_2" ] * coef[ "gamma_1" ]^( -coef[ "rho" ] ) * 
            ( exp( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
               ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) ) 
            )^( -coef[ "rho" ] ) +
            ( 1 - coef[ "delta_2" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] ) -
         ( coef[ "nu" ] / coef[ "rho" ] ) * 
         ( coef[ "delta_2" ] * coef[ "gamma_1" ]^( -coef[ "rho" ] ) * 
            ( exp( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
               ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) ) 
            )^( -coef[ "rho" ] ) +
            ( 1 - coef[ "delta_2" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
         ( - coef[ "delta_2" ] * log( coef[ "gamma_1" ] ) * 
            coef[ "gamma_1" ]^( -coef[ "rho" ] ) *
            ( exp( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
               ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) ) 
            )^( -coef[ "rho" ] ) -
            coef[ "delta_2" ] * coef[ "gamma_1" ]^( -coef[ "rho" ] ) *
            ( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
               ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) ) *
            ( exp( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
               ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) ) 
            )^( -coef[ "rho" ] ) -
            ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 3 ] ]] ) *
               data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) ) )
   } else {
      result <- -coef[ "gamma_2" ] * ( coef[ "nu" ] / coef[ "rho" ] ) * 
         B^( -coef[ "nu" ] / coef[ "rho" ] ) *
         ( -log( B ) / coef[ "rho" ] +
            B^(-1) * ( coef[ "delta_2" ] * coef[ "gamma_1" ]^( -coef[ "rho" ] ) * 
            B1^( coef[ "rho" ] / coef[ "rho_1" ] ) *
               ( - log( coef[ "gamma_1" ] ) + log( B1 ) / coef[ "rho_1" ] ) -
            ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 3 ] ]] ) *
            data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) ) )
   }

   return( result )
}


# derivatives with respect to nu
cesDerivCoefN3Nu <- function( coef, data, xNames ) {

   B1 <- cesDerivCoefN3B1( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN3B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         result <- coef[ "gamma_2" ] *
            ( coef[ "delta_2" ] *
               ( log( coef[ "gamma_1" ] ) + 
                  coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) ) +
               ( 1- coef[ "delta_2" ] ) * log( data[[ xNames[ 3 ] ]] ) ) *
            exp( - coef[ "nu" ] * ( - coef[ "delta_2" ] *
               ( log( coef[ "gamma_1" ] ) + 
                  coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) ) -
               ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      } else {
         result <- coef[ "gamma_2" ] * 
            ( coef[ "delta_2" ] *
               ( log( coef[ "gamma_1" ] ) - log( B1 ) / coef[ "rho_1" ] ) +
               ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 3 ] ]] ) ) *
            exp( - coef[ "nu" ] * ( - coef[ "delta_2" ] * 
                  ( log( coef[ "gamma_1" ] ) - log( B1 ) / coef[ "rho_1" ] ) -
               ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 3 ] ]] ) ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      result <- - ( coef[ "gamma_2" ] / coef[ "rho" ] ) *
         ( coef[ "delta_2" ] * coef[ "gamma_1" ]^( -coef[ "rho" ] ) *
            exp( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
               ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) 
            )^( -coef[ "rho" ] ) +
            ( 1 - coef[ "delta_2" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] )
         )^( - coef[ "nu" ] / coef[ "rho" ] ) *
         log( coef[ "delta_2" ] * coef[ "gamma_1" ]^( -coef[ "rho" ] ) * 
            exp( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) +
               ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) 
            )^( -coef[ "rho" ] ) +
            ( 1 - coef[ "delta_2" ] ) * data[[ xNames[ 3 ] ]]^( -coef[ "rho" ] ) )
   } else {
      result <- - coef[ "gamma_2" ] * log( B ) * 
         B^( -coef[ "nu" ] / coef[ "rho" ] ) / coef[ "rho" ]
   }

   return( result )
}

