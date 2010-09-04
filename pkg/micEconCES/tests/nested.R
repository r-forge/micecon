library( "micEconCES" )
library( maxLik )

data( "MishraCES" )

b <- c( "gamma" = 200, "delta_1" = 0.6, "delta_2" = 0.3, "rho_1" = 0.5,
   "rho_2" = -0.17, "rho" = 0.6 )
bVrs <- c( b, "nu" = 1.05 )
xNames <- c( "X1", "X2", "X3", "X4" )


## checking cesCalc()
MishraCES$Y2 <- cesCalc( xNames = xNames, 
   data = MishraCES, coef = b, nested = TRUE )

MishraCES$Y2

all.equal( MishraCES$Y, MishraCES$Y2 )

# VRS
MishraCES$yVrs <- cesCalc( xNames = xNames, 
   data = MishraCES, coef = bVrs, nested = TRUE )
MishraCES$yVrs


## checking cesDerivCoef
cesDeriv <- micEconCES:::cesDerivCoef( par = b, xNames = xNames, 
   data = MishraCES, vrs = FALSE, nested = TRUE )
f <- function( par ) {
   return( cesCalc( xNames = xNames, data = MishraCES, coef = par, 
      nested = TRUE ) )
}
cesDerivNum <- numericGradient( f, t0 = b )
all.equal( cesDeriv, cesDerivNum )

# VRE
cesDerivVrs <- micEconCES:::cesDerivCoef( par = bVrs, xNames = xNames, 
   data = MishraCES, vrs = TRUE, nested = TRUE )
cesDerivVrsNum <- numericGradient( f, t0 = bVrs )

all.equal( cesDerivVrs, cesDerivVrsNum )

