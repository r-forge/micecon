# calculate part "B1"
cesDerivCoefN4B1 <- function( par, data, xNames ) {

      B1 <- par[ "delta_1" ] * data[[ xNames[ 1 ] ]]^(-par[ "rho_1" ]) + 
         ( 1 - par[ "delta_1" ] ) * data[[ xNames[ 2 ] ]]^(-par[ "rho_1" ])

   return( B1 )
}


# calculate part "B2"
cesDerivCoefN4B2 <- function( par, data, xNames ) {

      B2 <- par[ "delta_2" ] * data[[ xNames[ 3 ] ]]^(-par[ "rho_2" ]) + 
         ( 1 - par[ "delta_2" ] ) * data[[ xNames[ 4 ] ]]^(-par[ "rho_2" ])

   return( B2 )
}


# calculate part "B"
cesDerivCoefN4B <- function( par, data, xNames ) {

   B1 <- cesDerivCoefN4B1( par = par, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( par = par, data = data, xNames = xNames )

      B <- par[ "delta_3" ] * B1^( par[ "rho" ] / par[ "rho_1" ] ) + 
         ( 1 - par[ "delta_3" ] ) * B2^( par[ "rho" ] / par[ "rho_2" ] )

   return( B )
}


# derivatives with respect to gamma
cesDerivCoefN4Gamma <- function( par, data, xNames ) {

   B <- cesDerivCoefN4B( par = par, data = data, xNames = xNames )

      result <- B^(-par[ "nu" ]/par[ "rho" ])

   return( result )
}


# derivatives with respect to delta_1
cesDerivCoefN4Delta1 <- function( par, data, xNames ) {

   B1 <- cesDerivCoefN4B1( par = par, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( par = par, data = data, xNames = xNames )

      result <- par[ "gamma" ] * 
         ( -par[ "nu" ] / par[ "rho" ] ) * 
         B^((-par[ "nu" ]-par[ "rho" ])/par[ "rho" ]) *
         (par[ "rho" ]/par[ "rho_1" ]) * par[ "delta_3" ] * 
         B1^((par[ "rho" ]-par[ "rho_1" ])/par[ "rho_1" ]) * 
         ( data[[ xNames[ 1 ] ]]^(-par[ "rho_1" ]) - 
            data[[ xNames[ 2 ] ]]^(-par[ "rho_1" ]) )

   return( result )
}


# derivatives with respect to delta_2
cesDerivCoefN4Delta2 <- function( par, data, xNames ) {

   B2 <- cesDerivCoefN4B2( par = par, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( par = par, data = data, xNames = xNames )

      result <- par[ "gamma" ] * 
         ( -par[ "nu" ] / par[ "rho" ] ) * 
         B^((-par[ "nu" ]-par[ "rho" ])/par[ "rho" ]) *
         (par[ "rho" ]/par[ "rho_2" ]) * ( 1 - par[ "delta_3" ] ) * 
         B2^((par[ "rho" ]-par[ "rho_2" ])/par[ "rho_2" ]) * 
         ( data[[ xNames[ 3 ] ]]^(-par[ "rho_2" ]) - 
            data[[ xNames[ 4 ] ]]^(-par[ "rho_2" ]) )

   return( result )
}

# derivatives with respect to delta_3
cesDerivCoefN4Delta3 <- function( par, data, xNames ) {

   B1 <- cesDerivCoefN4B1( par = par, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( par = par, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( par = par, data = data, xNames = xNames )

      result <- par[ "gamma" ] * 
         ( -par[ "nu" ] / par[ "rho" ] ) * 
         B^((-par[ "nu" ]-par[ "rho" ])/par[ "rho" ]) *
         ( B1^(par[ "rho" ]/par[ "rho_1" ]) - 
            B2^(par[ "rho" ]/par[ "rho_2" ]) )

   return( result )
}


# derivatives with respect to rho_1
cesDerivCoefN4Rho1 <- function( par, data, xNames ) {

   B1 <- cesDerivCoefN4B1( par = par, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( par = par, data = data, xNames = xNames )

         result <- par[ "gamma" ] * ( -par[ "nu" ] / par[ "rho" ] ) * 
            B^((-par[ "nu" ]-par[ "rho" ])/par[ "rho" ]) *
            ( par[ "delta_3" ] * log( B1 ) * B1^(par[ "rho" ]/par[ "rho_1" ]) * 
               ( -par[ "rho" ]/par[ "rho_1" ]^2 ) + 
               par[ "delta_3" ] * 
               B1^((par[ "rho" ]-par[ "rho_1" ])/par[ "rho_1" ]) * 
               (par[ "rho" ]/par[ "rho_1" ]) *
               ( -par[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) * 
                  data[[ xNames[ 1 ] ]]^(-par[ "rho_1" ]) - 
                  ( 1 - par[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) * 
                  data[[ xNames[ 2 ] ]]^(-par[ "rho_1" ]) ) )

   return( result )
}


# derivatives with respect to rho_2
cesDerivCoefN4Rho2 <- function( par, data, xNames ) {

   B2 <- cesDerivCoefN4B2( par = par, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( par = par, data = data, xNames = xNames )

         result <- par[ "gamma" ] * ( -par[ "nu" ] / par[ "rho" ] ) * 
            B^((-par[ "nu" ]-par[ "rho" ])/par[ "rho" ]) *
            ( ( 1 - par[ "delta_3" ] ) * log( B2 ) * 
               B2^(par[ "rho" ]/par[ "rho_2" ]) * 
               ( -par[ "rho" ]/par[ "rho_2" ]^2 ) + 
               ( 1 - par[ "delta_3" ] ) * 
               B2^((par[ "rho" ]-par[ "rho_2" ])/par[ "rho_2" ]) * 
               (par[ "rho" ]/par[ "rho_2" ]) *
               ( -par[ "delta_2" ] * log( data[[ xNames[ 3 ] ]] ) * 
                  data[[ xNames[ 3 ] ]]^(-par[ "rho_2" ]) - 
                  ( 1 - par[ "delta_2" ] ) * log( data[[ xNames[ 4 ] ]] ) * 
                  data[[ xNames[ 4 ] ]]^(-par[ "rho_2" ]) ) )

   return( result )
}

# derivatives with respect to rho
cesDerivCoefN4Rho <- function( par, data, xNames ) {

   B1 <- cesDerivCoefN4B1( par = par, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( par = par, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( par = par, data = data, xNames = xNames )

         result <- par[ "gamma" ] * log( B ) * 
            B^(-par[ "nu" ]/par[ "rho" ]) * ( par[ "nu" ] / par[ "rho" ]^2 ) +
            par[ "gamma" ] * ( -par[ "nu" ]/par[ "rho" ] ) * 
            B^((-par[ "nu" ]-par[ "rho" ])/par[ "rho" ]) * 
            ( par[ "delta_3" ] * log( B1 ) * 
               B1^(par[ "rho" ]/par[ "rho_1" ]) / par[ "rho_1" ] + 
               ( 1 - par[ "delta_3" ] ) * log( B2 ) * 
               B2^(par[ "rho" ]/par[ "rho_2" ]) / par[ "rho_2" ] )

   return( result )
}

# derivatives with respect to nu
cesDerivCoefN4Nu <- function( par, data, xNames ) {

   B <- cesDerivCoefN4B( par = par, data = data, xNames = xNames )

         result <- - par[ "gamma" ] * log( B ) * 
            B^( -par[ "nu" ] / par[ "rho" ] ) / par[ "rho" ]

   return( result )
}

