# calculate part "B1"
cesDerivCoefN4B1 <- function( coef, data, xNames ) {

      B1 <- coef[ "delta_1" ] * data[[ xNames[ 1 ] ]]^(-coef[ "rho_1" ]) + 
         ( 1 - coef[ "delta_1" ] ) * data[[ xNames[ 2 ] ]]^(-coef[ "rho_1" ])

   return( B1 )
}


# calculate part "BB1"
cesDerivCoefN4BB1 <- function( coef, data, xNames ) {

   BB1 <- coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) + 
      ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] )

   return( BB1 )
}


# calculate part "B2"
cesDerivCoefN4B2 <- function( coef, data, xNames ) {

      B2 <- coef[ "delta_2" ] * data[[ xNames[ 3 ] ]]^(-coef[ "rho_2" ]) + 
         ( 1 - coef[ "delta_2" ] ) * data[[ xNames[ 4 ] ]]^(-coef[ "rho_2" ])

   return( B2 )
}


# calculate part "BB2"
cesDerivCoefN4BB2 <- function( coef, data, xNames ) {

   BB2 <- coef[ "delta_2" ] * log( data[[ xNames[ 3 ] ]] ) + 
      ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 4 ] ]] )

   return( BB2 )
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

   BB1 <- cesDerivCoefN4BB1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   BB2 <- cesDerivCoefN4BB2( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         if( coef[ "rho_2" ] == 0 ) {
            result <- exp( - coef[ "nu" ] * ( coef[ "delta_3" ] * ( - BB1 ) +
                  ( 1 - coef[ "delta_3" ] ) * ( - BB2 ) ) )
         } else {
            result <- exp( - coef[ "nu" ] * ( coef[ "delta_3" ] * ( - BB1 ) +
                  ( 1- coef[ "delta_3" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
         }
      } else if( coef[ "rho_2" ] == 0 ) {
         result <- exp( - coef[ "nu" ] * ( coef[ "delta_3" ] *
               log( B1 ) / coef[ "rho_1" ] +
               ( 1 - coef[ "delta_3" ] ) * ( - BB2 ) ) )
      } else {
         result <- 
            exp( - coef[ "nu" ] *
               ( coef[ "delta_3" ] * log( B1 ) / coef[ "rho_1" ] +
                  ( 1 - coef[ "delta_3" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      if( coef[ "rho_2" ] == 0 ) {
         result <- exp( - ( coef[ "nu" ] / coef[ "rho" ] ) *
            log( coef[ "delta_3" ] * exp( coef[ "rho" ] * ( - BB1 ) ) +
            ( 1 - coef[ "delta_3" ] ) * exp( coef[ "rho" ] * ( - BB2 ) )
            ) )
      } else {
         result <- 
            ( coef[ "delta_3" ] * 
               exp( coef[ "rho" ] * ( - BB1 ) ) +
               ( 1 - coef[ "delta_3" ] ) * 
                  B2^( coef[ "rho" ] / coef[ "rho_2" ] ) 
            )^( - coef[ "nu" ] / coef[ "rho" ] )
      }
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

   BB1 <- cesDerivCoefN4BB1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   BB2 <- cesDerivCoefN4BB2( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         result <- coef[ "gamma" ] *
            ( - coef[ "nu" ]* 
            ( coef[ "delta_3" ] *
               ( - log( data[[ xNames[ 1 ] ]] ) +
                  log( data[[ xNames[ 2 ] ]] ) ) ) ) *
            exp( -coef[ "nu" ] * 
               ( coef[ "delta_3" ] * ( - BB1 ) +
               ( 1 - coef[ "delta_3" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
      } else if( coef[ "rho_2" ] == 0 ) {
         result <- coef[ "gamma" ] * 
            ( -coef[ "nu" ] * coef[ "delta_3" ] * 
               ( data[[ xNames[ 1 ] ]]^( - coef[ "rho_1" ] ) -
                  data[[ xNames[ 2 ] ]]^( - coef[ "rho_1" ] ) ) ) /
            ( coef[ "rho_1" ] * B1 ) *
            exp( -coef[ "nu" ] *
               ( coef[ "delta_3" ] * log( B1 ) / coef[ "rho_1" ] +
                  ( 1 - coef[ "delta_3" ] ) * ( - BB2 ) ) )
      } else {
         result <- coef[ "gamma" ] *
            ( - coef[ "nu" ] / coef[ "rho_1" ] ) * 
            ( coef[ "delta_3" ] *
               ( data[[ xNames[ 1 ] ]]^(-coef[ "rho_1" ]) - 
                  data[[ xNames[ 2 ] ]]^(-coef[ "rho_1" ]) ) / B1 ) *
            exp( -coef[ "nu" ] * 
               ( coef[ "delta_3" ] * log( B1 ) / coef[ "rho_1" ] +
                  ( 1 - coef[ "delta_3" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      if( coef[ "rho_2" ] == 0 ) {
         result <- coef[ "gamma" ] * ( -coef[ "nu" ] / coef[ "rho" ] ) *
            ( coef[ "delta_3" ] * exp( coef[ "rho" ] * ( - BB1 ) ) +
               ( 1 - coef[ "delta_3" ] ) * exp( coef[ "rho" ] * ( - BB2 ) )
            )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
            ( coef[ "delta_3" ] * exp( coef[ "rho" ] * ( - BB1 ) ) *
                  coef[ "rho" ] * ( - log( data[[ xNames[ 1 ] ]] ) + 
                     log( data[[ xNames[ 2 ] ]] ) ) )
      } else {
         result <- coef[ "gamma" ] * 
            ( -coef[ "nu" ] / coef[ "rho" ] ) *
            ( coef[ "delta_3" ] *
               exp( coef[ "rho" ] * ( - BB1 ) ) +
               ( 1 - coef[ "delta_3" ] ) * 
                  B2^( coef[ "rho" ] / coef[ "rho_2" ] ) 
            )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
            ( coef[ "delta_3" ] *
               exp( coef[ "rho" ] *
                  ( - BB1 ) ) *
            coef[ "rho" ] *
            ( - log( data[[ xNames[ 1 ] ]] ) + log( data[[ xNames[ 2 ] ]] ) ) )
      }
   } else if( coef[ "rho_2" ] == 0 ) {
      result <- coef[ "gamma" ] * 
         ( -coef[ "nu" ] / coef[ "rho" ] ) *
         ( ( 1 - coef[ "delta_3" ] ) *
            exp( coef[ "rho" ] * ( - BB2 ) ) +
            coef[ "delta_3" ] * 
               B1^( coef[ "rho" ] / coef[ "rho_1" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
         ( coef[ "delta_3" ] * coef[ "rho" ] *
            B1^( coef[ "rho" ] / coef[ "rho_1" ] - 1 ) *
            ( data[[ xNames[ 1 ] ]]^(-coef[ "rho_1" ]) - 
               data[[ xNames[ 2 ] ]]^(-coef[ "rho_1" ]) ) /
            coef[ "rho_1" ] )
   } else {
      result <- coef[ "gamma" ] * 
         ( -coef[ "nu" ] / coef[ "rho" ] ) * 
         B^((-coef[ "nu" ]-coef[ "rho" ])/coef[ "rho" ]) *
         (coef[ "rho" ]/coef[ "rho_1" ]) * coef[ "delta_3" ] * 
         B1^((coef[ "rho" ]-coef[ "rho_1" ])/coef[ "rho_1" ]) * 
         ( data[[ xNames[ 1 ] ]]^(-coef[ "rho_1" ]) - 
            data[[ xNames[ 2 ] ]]^(-coef[ "rho_1" ]) )
   }

   return( result )
}


# derivatives with respect to delta_2
cesDerivCoefN4Delta2 <- function( coef, data, xNames ) {

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   BB1 <- cesDerivCoefN4BB1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   BB2 <- cesDerivCoefN4BB2( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         result <- coef[ "gamma" ] * 
            ( -coef[ "nu" ] * ( 1 - coef[ "delta_3" ] ) * 
               ( data[[ xNames[ 3 ] ]]^( - coef[ "rho_2" ] ) -
                  data[[ xNames[ 4 ] ]]^( - coef[ "rho_2" ] ) ) ) /
            ( coef[ "rho_2" ] * B2 ) *
            exp( -coef[ "nu" ] *
               ( coef[ "delta_3" ] * ( -BB1 ) +
                  ( 1 - coef[ "delta_3" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
      } else if( coef[ "rho_2" ] == 0 ) {
         result <- coef[ "gamma" ] *
            ( - coef[ "nu" ]* 
            ( ( 1 - coef[ "delta_3" ] ) *
               ( - log( data[[ xNames[ 3 ] ]] ) +
                  log( data[[ xNames[ 4 ] ]] ) ) ) ) *
            exp( -coef[ "nu" ] * 
               ( coef[ "delta_3" ] * log( B1 ) / coef[ "rho_1" ] +
               ( 1 - coef[ "delta_3" ] ) * ( - BB2 ) ) )
      } else {
         result <- coef[ "gamma" ] *
            ( - coef[ "nu" ] / coef[ "rho_2" ] ) * 
            ( ( 1 - coef[ "delta_3" ] ) *
               ( data[[ xNames[ 3 ] ]]^(-coef[ "rho_2" ]) - 
                  data[[ xNames[ 4 ] ]]^(-coef[ "rho_2" ]) ) / B2 ) *
            exp( -coef[ "nu" ] * 
               ( coef[ "delta_3" ] * log( B1 ) / coef[ "rho_1" ] +
                  ( 1 - coef[ "delta_3" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      if( coef[ "rho_2" ] == 0 ) {
         result <- - coef[ "gamma" ] * ( coef[ "nu" ] / coef[ "rho" ] ) *
            ( coef[ "delta_3" ] * exp( coef[ "rho" ] * ( - BB1 ) ) +
               ( 1 - coef[ "delta_3" ] ) * exp( coef[ "rho" ] * ( - BB2 ) )
            )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
            ( ( 1 - coef[ "delta_3" ] ) * exp( coef[ "rho" ] * ( - BB2 ) ) *
                  coef[ "rho" ] * ( - log( data[[ xNames[ 3 ] ]] ) + 
                     log( data[[ xNames[ 4 ] ]] ) ) )
      } else {
         result <- coef[ "gamma" ] * 
            ( -coef[ "nu" ] / coef[ "rho" ] ) *
            ( coef[ "delta_3" ] *
               exp( coef[ "rho" ] * ( - BB1 ) ) +
               ( 1 - coef[ "delta_3" ] ) * 
                  B2^( coef[ "rho" ] / coef[ "rho_2" ] ) 
            )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
            ( ( 1 - coef[ "delta_3" ] ) * coef[ "rho" ] *
               B2^( coef[ "rho" ] / coef[ "rho_2" ] - 1 ) *
               ( data[[ xNames[ 3 ] ]]^(-coef[ "rho_2" ]) - 
                  data[[ xNames[ 4 ] ]]^(-coef[ "rho_2" ]) ) /
               coef[ "rho_2" ] )
      }
   } else if( coef[ "rho_2" ] == 0 ) {
      result <- coef[ "gamma" ] * 
         ( -coef[ "nu" ] / coef[ "rho" ] ) *
         ( ( 1 - coef[ "delta_3" ] ) *
            exp( coef[ "rho" ] * ( - BB2 ) ) +
            coef[ "delta_3" ] * 
               B1^( coef[ "rho" ] / coef[ "rho_1" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
         ( ( 1 - coef[ "delta_3" ] ) *
            exp( coef[ "rho" ] * ( - BB2 ) ) *
         coef[ "rho" ] *
         ( - log( data[[ xNames[ 3 ] ]] ) + log( data[[ xNames[ 4 ] ]] ) ) )
   } else {
      result <- coef[ "gamma" ] * 
         ( -coef[ "nu" ] / coef[ "rho" ] ) * 
         B^((-coef[ "nu" ]-coef[ "rho" ])/coef[ "rho" ]) *
         (coef[ "rho" ]/coef[ "rho_2" ]) * ( 1 - coef[ "delta_3" ] ) * 
         B2^((coef[ "rho" ]-coef[ "rho_2" ])/coef[ "rho_2" ]) * 
         ( data[[ xNames[ 3 ] ]]^(-coef[ "rho_2" ]) - 
            data[[ xNames[ 4 ] ]]^(-coef[ "rho_2" ]) )
   }

   return( result )
}

# derivatives with respect to delta_3
cesDerivCoefN4Delta3 <- function( coef, data, xNames ) {

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   BB1 <- cesDerivCoefN4BB1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   BB2 <- cesDerivCoefN4BB2( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      result <- coef[ "gamma" ] *
         ( - coef[ "nu" ] * 
            ( log( B1 ) / coef[ "rho_1" ] - log( B2 ) / coef[ "rho_2" ] ) ) *
         exp( -coef[ "nu" ] * 
            ( coef[ "delta_3" ] * log( B1 ) / coef[ "rho_1" ] +
               ( 1 - coef[ "delta_3" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
   } else if( coef[ "rho_1" ] == 0 ) {
      if( coef[ "rho_2" ] == 0 ) {
         result <- coef[ "gamma" ] * 
            ( -coef[ "nu" ] / coef[ "rho" ] ) *
            ( coef[ "delta_3" ] * exp( coef[ "rho" ] * ( - BB1 ) ) +
               ( 1 - coef[ "delta_3" ] ) * exp( coef[ "rho" ] * ( - BB2 ) ) 
            )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
            ( exp( coef[ "rho" ] * ( - BB1 ) ) -
               exp( coef[ "rho" ] * ( - BB2 ) )
            ) 
      } else {
         result <- coef[ "gamma" ] * 
            ( -coef[ "nu" ] / coef[ "rho" ] ) *
            ( coef[ "delta_3" ] *
               exp( coef[ "rho" ] * ( - BB1 ) ) +
               ( 1 - coef[ "delta_3" ] ) * 
                  B2^( coef[ "rho" ] / coef[ "rho_2" ] ) 
            )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
            ( exp( coef[ "rho" ] * ( - BB1 ) ) -
               B2^( coef[ "rho" ] / coef[ "rho_2" ] ) )
      }
   } else if( coef[ "rho_2" ] == 0 ) {
      result <- coef[ "gamma" ] * 
         ( -coef[ "nu" ] / coef[ "rho" ] ) *
         ( ( 1 -coef[ "delta_3" ] ) *
            exp( coef[ "rho" ] * ( - BB2 ) ) +
            coef[ "delta_3" ] * 
               B1^( coef[ "rho" ] / coef[ "rho_1" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
         ( - exp( coef[ "rho" ] * ( - BB2 ) ) +
            B1^( coef[ "rho" ] / coef[ "rho_1" ] ) )
   } else {
      result <- coef[ "gamma" ] * 
         ( -coef[ "nu" ] / coef[ "rho" ] ) * 
         B^((-coef[ "nu" ]-coef[ "rho" ])/coef[ "rho" ]) *
         ( B1^(coef[ "rho" ]/coef[ "rho_1" ]) - 
            B2^(coef[ "rho" ]/coef[ "rho_2" ]) )
   }

   return( result )
}


# derivatives with respect to rho_1
cesDerivCoefN4Rho1 <- function( coef, data, xNames ) {

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   BB1 <- cesDerivCoefN4BB1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   BB2 <- cesDerivCoefN4BB2( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      result <- coef[ "gamma" ] * 
         ( - coef[ "nu" ] *
            ( coef[ "delta_3" ] * coef[ "rho_1" ] *
               ( -coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) * 
                  data[[ xNames[ 1 ] ]]^(-coef[ "rho_1" ] ) - 
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) * 
                  data[[ xNames[ 2 ] ]]^(-coef[ "rho_1" ] ) ) /
                  ( B1 * coef[ "rho_1" ]^2   ) -
                  coef[ "delta_3" ] * log( B1 ) / coef[ "rho_1" ]^2 ) ) *
         exp( - coef[ "nu" ] *
            ( coef[ "delta_3" ] * log( B1 ) / coef[ "rho_1" ] +
               ( 1 - coef[ "delta_3" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
   } else if( coef[ "rho_1" ] == 0 ) {
      if( coef[ "rho_2" ] == 0 ) {
         result <- - coef[ "gamma" ] * coef[ "delta_3" ] * coef[ "nu" ] * 
            ( coef[ "delta_3" ] *  exp( coef[ "rho" ] * ( - BB1 ) ) +
               ( 1 - coef[ "delta_3" ] ) * exp( coef[ "rho" ] * ( - BB2 ) )
            )^( -coef[ "nu" ] / coef[ "rho" ] - 1 ) *
            exp( - coef[ "rho" ] * BB1 ) *
                  ( - 0.5 * BB1^2 +
                     0.5 * ( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] )^2 +
                        ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] )^2 )
                  )
      } else {
         result <- - coef[ "gamma" ] * coef[ "nu" ] * coef[ "delta_3" ] * 
            ( coef[ "delta_3" ] *
               exp( coef[ "rho" ] * ( - BB1 ) ) +
               ( 1 - coef[ "delta_3" ] ) * 
               B2^( coef[ "rho" ] / coef[ "rho_2" ] )
            )^( -coef[ "nu" ] / coef[ "rho" ] - 1 ) *
            exp( - coef[ "rho" ] * BB1 ) *
            ( - 0.5 * BB1^2 +
               0.5 * ( coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] )^2 +
                  ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] )^2 )
            )
      }
   } else if( coef[ "rho_2" ] == 0 ) {
      result <- coef[ "gamma" ] * ( -coef[ "nu" ] / coef[ "rho" ] ) * 
         ( ( 1 - coef[ "delta_3" ] ) *
            exp( coef[ "rho" ] * ( - BB2 ) ) +
            coef[ "delta_3" ] * 
            B1^( coef[ "rho" ] / coef[ "rho_1" ] )
         )^( -coef[ "nu" ] / coef[ "rho" ] - 1 ) *
         ( coef[ "delta_3" ] * log( B1 ) * 
            B1^( coef[ "rho" ] / coef[ "rho_1" ] ) *
            ( - coef[ "rho" ] / coef[ "rho_1" ]^2 ) +
            coef[ "delta_3" ] * ( coef[ "rho" ] / coef[ "rho_1" ] ) *
            B1^( coef[ "rho" ] / coef[ "rho_1" ] - 1 ) *
            ( - coef[ "delta_1" ] * log( data[[ xNames[ 1 ] ]] ) *
               data[[ xNames[ 1 ] ]]^( - coef[ "rho_1" ] ) -
               ( 1 - coef[ "delta_1" ] ) * log( data[[ xNames[ 2 ] ]] ) *
               data[[ xNames[ 2 ] ]]^( - coef[ "rho_1" ] )
            )
         )
   } else {
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
   }

   return( result )
}


# derivatives with respect to rho_2
cesDerivCoefN4Rho2 <- function( coef, data, xNames ) {

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   BB1 <- cesDerivCoefN4BB1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   BB2 <- cesDerivCoefN4BB2( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      result <- coef[ "gamma" ] * 
         ( - coef[ "nu" ] *
            ( ( 1 - coef[ "delta_3" ] ) * coef[ "rho_2" ] *
               ( -coef[ "delta_2" ] * log( data[[ xNames[ 3 ] ]] ) * 
                  data[[ xNames[ 3 ] ]]^(-coef[ "rho_2" ] ) - 
                  ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 4 ] ]] ) * 
                  data[[ xNames[ 4 ] ]]^(-coef[ "rho_2" ] ) ) /
                  ( B2 * coef[ "rho_2" ]^2   ) -
                  ( 1 - coef[ "delta_3" ] ) * log( B2 ) / coef[ "rho_2" ]^2 ) ) *
         exp( - coef[ "nu" ] *
            ( coef[ "delta_3" ] * log( B1 ) / coef[ "rho_1" ] +
               ( 1 - coef[ "delta_3" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
   } else if( coef[ "rho_1" ] == 0 ) {
      result <- coef[ "gamma" ] * ( -coef[ "nu" ] / coef[ "rho" ] ) *
         ( coef[ "delta_3" ] *
            exp( coef[ "rho" ] * ( - BB1 ) ) +
            ( 1 - coef[ "delta_3" ] ) * 
               B2^( coef[ "rho" ] / coef[ "rho_2" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
         ( ( 1 - coef[ "delta_3" ] ) * log( B2 ) * 
            B2^( coef[ "rho" ] / coef[ "rho_2" ] ) *
            ( - coef[ "rho" ] / coef[ "rho_2" ]^2 ) +
            ( 1 - coef[ "delta_3" ] ) * ( coef[ "rho" ] / coef[ "rho_2" ] ) *
            B2^( coef[ "rho" ] / coef[ "rho_2" ] - 1 ) *
            ( -coef[ "delta_2" ] * log( data[[ xNames[ 3 ] ]] ) * 
               data[[ xNames[ 3 ] ]]^( - coef[ "rho_2" ]) - 
               ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 4 ] ]] ) * 
               data[[ xNames[ 4 ] ]]^( - coef[ "rho_2" ] ) ) )
   } else if( coef[ "rho_2" ] == 0 ) {
      result <- - coef[ "gamma" ] * coef[ "nu" ] * ( 1 - coef[ "delta_3" ] ) * 
         ( ( 1 - coef[ "delta_3" ] ) *
            exp( coef[ "rho" ] * ( - BB2 ) ) +
            coef[ "delta_3" ] * 
            B1^( coef[ "rho" ] / coef[ "rho_1" ] )
         )^( -coef[ "nu" ] / coef[ "rho" ] - 1 ) *
         exp( - coef[ "rho" ] * BB2 ) *
         ( - 0.5 * BB2^2 +
            0.5 * ( coef[ "delta_2" ] * log( data[[ xNames[ 3 ] ]] )^2 +
               ( 1 - coef[ "delta_2" ] ) * log( data[[ xNames[ 4 ] ]] )^2 )
         )
   } else {
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
   }

   return( result )
}

# derivatives with respect to rho
cesDerivCoefN4Rho <- function( coef, data, xNames ) {

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   BB1 <- cesDerivCoefN4BB1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   BB2 <- cesDerivCoefN4BB2( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      result <- coef[ "gamma" ] * coef[ "nu" ] *
         ( -0.5 * ( coef[ "delta_3" ] * log( B1 )^2 / coef[ "rho_1" ]^2 +
               ( 1 - coef[ "delta_3" ] ) * log( B2 )^2 / coef[ "rho_2" ]^2 ) +
            0.5 * ( coef[ "delta_3" ] * log( B1 ) / coef[ "rho_1" ] +
               ( 1 - coef[ "delta_3" ] ) * log( B2 ) / coef[ "rho_2" ] )^2 ) *
         exp( - coef[ "nu" ] * 
            ( coef[ "delta_3" ] * log( B1 ) / coef[ "rho_1" ] +
               ( 1 - coef[ "delta_3" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
   } else if( coef[ "rho_1" ] == 0 ) {
      result <- coef[ "gamma" ] *
         log( coef[ "delta_3" ] *
            exp( coef[ "rho" ] * ( - BB1 ) ) +
            ( 1 - coef[ "delta_3" ] ) * 
               B2^( coef[ "rho" ] / coef[ "rho_2" ] ) ) *
         ( coef[ "delta_3" ] *
            exp( coef[ "rho" ] * ( - BB1 ) ) +
            ( 1 - coef[ "delta_3" ] ) * 
               B2^( coef[ "rho" ] / coef[ "rho_2" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] ) *
         ( coef[ "nu" ] / coef[ "rho" ]^2 ) +
         coef[ "gamma" ] * ( -coef[ "nu" ] / coef[ "rho" ] ) *
         ( coef[ "delta_3" ] *
            exp( coef[ "rho" ] * ( - BB1 ) ) +
            ( 1 - coef[ "delta_3" ] ) * 
               B2^( coef[ "rho" ] / coef[ "rho_2" ] )
         )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
         ( coef[ "delta_3" ] *
            exp( coef[ "rho" ] * ( - BB1 ) ) * ( - BB1 ) +
            ( 1 - coef[ "delta_3" ] ) * log( B2 ) * 
               B2^( coef[ "rho" ] / coef[ "rho_2" ] ) / coef[ "rho_2" ] )
   } else if( coef[ "rho_2" ] == 0 ) {
      result <- coef[ "gamma" ] *
         log( ( 1 - coef[ "delta_3" ] ) *
            exp( coef[ "rho" ] * ( - BB2 ) ) +
            coef[ "delta_3" ] *
               B1^( coef[ "rho" ] / coef[ "rho_1" ] ) ) *
         ( ( 1 - coef[ "delta_3" ] ) *
            exp( coef[ "rho" ] * ( - BB2 ) ) +
            coef[ "delta_3" ] * 
               B1^( coef[ "rho" ] / coef[ "rho_1" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] ) *
         ( coef[ "nu" ] / coef[ "rho" ]^2 ) +
         coef[ "gamma" ] * ( -coef[ "nu" ] / coef[ "rho" ] ) *
         ( ( 1 - coef[ "delta_3" ] ) *
            exp( coef[ "rho" ] * ( - BB2 ) ) +
            coef[ "delta_3" ] *
               B1^( coef[ "rho" ] / coef[ "rho_1" ] )
         )^( - coef[ "nu" ] / coef[ "rho" ] - 1 ) *
         ( ( 1 - coef[ "delta_3" ] ) *
            exp( coef[ "rho" ] * ( - BB2 ) ) * ( - BB2 ) +
            ( coef[ "delta_3" ] * log( B1 ) * 
               B1^( coef[ "rho" ] / coef[ "rho_1" ] ) / coef[ "rho_1" ] )
         )
   } else {
         result <- coef[ "gamma" ] * log( B ) * 
            B^(-coef[ "nu" ]/coef[ "rho" ]) * ( coef[ "nu" ] / coef[ "rho" ]^2 ) +
            coef[ "gamma" ] * ( -coef[ "nu" ]/coef[ "rho" ] ) * 
            B^((-coef[ "nu" ]-coef[ "rho" ])/coef[ "rho" ]) * 
            ( coef[ "delta_3" ] * log( B1 ) * 
               B1^(coef[ "rho" ]/coef[ "rho_1" ]) / coef[ "rho_1" ] + 
               ( 1 - coef[ "delta_3" ] ) * log( B2 ) * 
               B2^(coef[ "rho" ]/coef[ "rho_2" ]) / coef[ "rho_2" ] )
   }

   return( result )
}

# derivatives with respect to nu
cesDerivCoefN4Nu <- function( coef, data, xNames ) {

   B1 <- cesDerivCoefN4B1( coef = coef, data = data, xNames = xNames )

   BB1 <- cesDerivCoefN4BB1( coef = coef, data = data, xNames = xNames )

   B2 <- cesDerivCoefN4B2( coef = coef, data = data, xNames = xNames )

   BB2 <- cesDerivCoefN4BB2( coef = coef, data = data, xNames = xNames )

   B <- cesDerivCoefN4B( coef = coef, data = data, xNames = xNames )

   if( coef[ "rho" ] == 0 ) {
      if( coef[ "rho_1" ] == 0 ) {
         result <- coef[ "gamma" ] *
            ( - coef[ "delta_3" ] * ( - BB1 ) -
            ( 1 - coef[ "delta_3" ] ) * log( B2 ) / coef[ "rho_2" ] ) *
            exp( - coef[ "nu" ] *  
               ( coef[ "delta_3" ] * ( - BB1 ) +
                  ( 1 - coef[ "delta_3" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
      } else {
         result <- - coef[ "gamma" ] * 
            ( coef[ "delta_3" ] * log( B1 ) / coef[ "rho_1" ] +
               ( 1 - coef[ "delta_3" ] ) * log( B2 ) / coef[ "rho_2" ] ) *
            exp( - coef[ "nu" ] *
               ( coef[ "delta_3" ] * log( B1 ) / coef[ "rho_1" ] +
                  ( 1 - coef[ "delta_3" ] ) * log( B2 ) / coef[ "rho_2" ] ) )
      }
   } else if( coef[ "rho_1" ] == 0 ) {
      if( coef[ "rho_2" ] == 0 ) {
         result <- - ( coef[ "gamma" ] / coef[ "rho" ] ) * 
            log( ( 1 - coef[ "delta_3" ] ) * exp( coef[ "rho" ] * ( - BB2 ) ) +
               coef[ "delta_3" ] * exp( coef[ "rho" ] * ( - BB1 ) )
            ) *            
            ( coef[ "delta_3" ] * exp( coef[ "rho" ] * ( - BB1 ) ) +
               ( 1 - coef[ "delta_3" ] ) * exp( coef[ "rho" ] * ( - BB2 ) )
            )^( - coef[ "nu" ] / coef[ "rho" ] )
      } else {
         result <- - ( coef[ "gamma" ] / coef[ "rho" ] ) * 
            log( coef[ "delta_3" ] * 
               exp( coef[ "rho" ] * ( - BB1 ) ) +
               ( 1 - coef[ "delta_3" ] ) * 
               B2^( coef[ "rho" ] / coef[ "rho_2" ] ) ) *
            ( coef[ "delta_3" ] * 
               exp( coef[ "rho" ] * ( - BB1 ) ) +
               ( 1 - coef[ "delta_3" ] ) * 
               B2^( coef[ "rho" ] / coef[ "rho_2" ] ) 
            )^( - coef[ "nu" ] / coef[ "rho" ] )
      }
   } else if( coef[ "rho_2" ] == 0 ) {
      result <- - ( coef[ "gamma" ] / coef[ "rho" ] ) * 
         log( ( 1 - coef[ "delta_3" ] ) * 
            exp( coef[ "rho" ] * ( - BB2 ) ) +
            coef[ "delta_3" ] * 
            B1^( coef[ "rho" ] / coef[ "rho_1" ] ) ) *
         ( ( 1 - coef[ "delta_3" ] ) * 
            exp( coef[ "rho" ] * ( - BB2 ) ) +
            coef[ "delta_3" ] * 
            B1^( coef[ "rho" ] / coef[ "rho_1" ] ) 
         )^( - coef[ "nu" ] / coef[ "rho" ] )
   } else {
         result <- - coef[ "gamma" ] * log( B ) * 
            B^( -coef[ "nu" ] / coef[ "rho" ] ) / coef[ "rho" ]
   }

   return( result )
}

