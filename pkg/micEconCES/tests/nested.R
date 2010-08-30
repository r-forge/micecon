library( "micEconCES" )

data( "MishraCES" )

b <- c( "gamma" = 200, "delta_1" = 0.6, "delta_2" = 0.3, "rho_1" = 0.5,
   "rho_2" = -0.17, "rho" = 0.6 )

MishraCES$Y2 <- cesCalc( xNames = c( "X1", "X2", "X3", "X4" ), 
   data = MishraCES, coef = b, nested = TRUE )

MishraCES$Y2

all.equal( MishraCES$Y, MishraCES$Y2 )
