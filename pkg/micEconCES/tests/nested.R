library( "micEconCES" )

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
