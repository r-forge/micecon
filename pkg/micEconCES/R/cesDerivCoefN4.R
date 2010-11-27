# calculate part "B1"
cesDerivCoefN4B1 <- function( coef, data, xNames ) {

      B1 <- coef[ "delta_1" ] * data[[ xNames[ 1 ] ]]^(-coef[ "rho_1" ]) + 
         ( 1 - coef[ "delta_1" ] ) * data[[ xNames[ 2 ] ]]^(-coef[ "rho_1" ])

   return( B1 )
}


# calculate part "B2"
cesDerivCoefN4B2 <- function( coef, data, xNames ) {

      B2 <- coef[ "delta_2" ] * data[[ xNames[ 3 ] ]]^(-coef[ "rho_2" ]) + 
         ( 1 - coef[ "delta_2" ] ) * data[[ xNames[ 4 ] ]]^(-coef[ "rho_2" ])

   return( B2 )
}


# calculate part "B"
cesDerivCoefN4B <- function( coef, data, xNames ) {

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

      B <- coef[ "delta_3" ] * B1^( coef[ "rho" ] / coef[ "rho_1" ] ) + 
         ( 1 - coef[ "delta_3" ] ) * B2^( coef[ "rho" ] / coef[ "rho_2" ] )

   return( B )
}


# derivatives with respect to gamma
cesDerivCoefN4Gamma <- function( coef, data, xNames ) {

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho_1" ] == 0 ) {
      result <- 
         ( coef[ "delta_3" ] * 
            exp( coef[ "rho" ] * 
               ( - coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) -
               ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) ) ) +
            ( 1 - coef[ "delta_3" ] ) * 
               B2^( coef[ "rho" ] / coef[ "rho_2" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] )
   } else if( coef[ "rho_2" ] == 0 ) {
      result <-
         ( coef[ "delta_3" ] * B1^( coef[ "rho" ] / coef[ "rho_1" ] ) +
            ( 1 - coef[ "delta_3" ] ) * 
               exp( coef[ "rho" ] * 
                  ( - coef[ "delta_2" ] * log( data[[ xNames[ 3 ] ]] ) -
                     ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 4 ] ]] ) ) ) 
         )^( - ( coef[ "nu" ] / coef[ "rho" ] ) )
   } else {
      result <- B^(-coef[ "nu" ]/coef[ "rho" ])
   }

   return( result )
}


# derivatives with respect to delta_1
cesDerivCoefN4Delta1 <- function( coef, data, xNames ) {

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

      result <- coef[ "gamma" ] * 
         ( -coef[ "nu" ] / coef[ "rho" ] ) * 
         B^((-coef[ "nu" ]-coef[ "rho" ])/coef[ "rho" ]) *
         (coef[ "rho" ]/coef[ "rho_1" ]) * coef[ "delta_3" ] * 
         B1^((coef[ "rho" ]-coef[ "rho_1" ])/coef[ "rho_1" ]) * 
         ( data[[ xNames[ 1 ] ]]^(-coef[ "rho_1" ]) - 
            data[[ xNames[ 2 ] ]]^(-coef[ "rho_1" ]) )

   return( result )
}


# derivatives with respect to delta_2
cesDerivCoefN4Delta2 <- function( coef, data, xNames ) {

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

      result <- coef[ "gamma" ] * 
         ( -coef[ "nu" ] / coef[ "rho" ] ) * 
         B^((-coef[ "nu" ]-coef[ "rho" ])/coef[ "rho" ]) *
         (coef[ "rho" ]/coef[ "rho_2" ]) * ( 1 - coef[ "delta_3" ] ) * 
         B2^((coef[ "rho" ]-coef[ "rho_2" ])/coef[ "rho_2" ]) * 
         ( data[[ xNames[ 3 ] ]]^(-coef[ "rho_2" ]) - 
            data[[ xNames[ 4 ] ]]^(-coef[ "rho_2" ]) )

   return( result )
}

# derivatives with respect to delta_3
cesDerivCoefN4Delta3 <- function( coef, data, xNames ) {

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

      result <- coef[ "gamma" ] * 
         ( -coef[ "nu" ] / coef[ "rho" ] ) * 
         B^((-coef[ "nu" ]-coef[ "rho" ])/coef[ "rho" ]) *
         ( B1^(coef[ "rho" ]/coef[ "rho_1" ]) - 
            B2^(coef[ "rho" ]/coef[ "rho_2" ]) )

   return( result )
}


# derivatives with respect to rho_1
cesDerivCoefN4Rho1 <- function( coef, data, xNames ) {

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

         result <- coef[ "gamma" ] * ( -coef[ "nu" ] / coef[ "rho" ] ) * 
            B^((-coef[ "nu" ]-coef[ "rho" ])/coef[ "rho" ]) *
            ( coef[ "delta_3" ] * log( B1 ) * B1^(coef[ "rho" ]/coef[ "rho_1" ]) * 
               ( -coef[ "rho" ]/coef[ "rho_1" ]^2 ) + 
               coef[ "delta_3" ] * 
               B1^((coef[ "rho" ]-coef[ "rho_1" ])/coef[ "rho_1" ]) * 
               (coef[ "rho" ]/coef[ "rho_1" ]) *
               ( -coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) * 
                  data[[ xNames[ 1 ] ]]^(-coef[ "rho_1" ]) - 
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) * 
                  data[[ xNames[ 2 ] ]]^(-coef[ "rho_1" ]) ) )

   return( result )
}


# derivatives with respect to rho_2
cesDerivCoefN4Rho2 <- function( coef, data, xNames ) {

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

         result <- coef[ "gamma" ] * ( -coef[ "nu" ] / coef[ "rho" ] ) * 
            B^((-coef[ "nu" ]-coef[ "rho" ])/coef[ "rho" ]) *
            ( ( 1 - coef[ "delta_3" ] ) * log( B2 ) * 
               B2^(coef[ "rho" ]/coef[ "rho_2" ]) * 
               ( -coef[ "rho" ]/coef[ "rho_2" ]^2 ) + 
               ( 1 - coef[ "delta_3" ] ) * 
               B2^((coef[ "rho" ]-coef[ "rho_2" ])/coef[ "rho_2" ]) * 
               (coef[ "rho" ]/coef[ "rho_2" ]) *
               ( -coef[ "delta_2" ] * log( data[[ xNames[ 3 ] ]] ) * 
                  data[[ xNames[ 3 ] ]]^(-coef[ "rho_2" ]) - 
                  ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 4 ] ]] ) * 
                  data[[ xNames[ 4 ] ]]^(-coef[ "rho_2" ]) ) )

   return( result )
}

# derivatives with respect to rho
cesDerivCoefN4Rho <- function( coef, data, xNames ) {

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

         result <- coef[ "gamma" ] * log( B ) * 
            B^(-coef[ "nu" ]/coef[ "rho" ]) * ( coef[ "nu" ] / coef[ "rho" ]^2 ) +
            coef[ "gamma" ] * ( -coef[ "nu" ]/coef[ "rho" ] ) * 
            B^((-coef[ "nu" ]-coef[ "rho" ])/coef[ "rho" ]) * 
            ( coef[ "delta_3" ] * log( B1 ) * 
               B1^(coef[ "rho" ]/coef[ "rho_1" ]) / coef[ "rho_1" ] + 
               ( 1 - coef[ "delta_3" ] ) * log( B2 ) * 
               B2^(coef[ "rho" ]/coef[ "rho_2" ]) / coef[ "rho_2" ] )

   return( result )
}

# derivatives with respect to nu
cesDerivCoefN4Nu <- function( coef, data, xNames ) {

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

         result <- - coef[ "gamma" ] * log( B ) * 
            B^( -coef[ "nu" ] / coef[ "rho" ] ) / coef[ "rho" ]

   return( result )
}

