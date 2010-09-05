library( "micEconCES" )

data( "MishraCES" )

b <- c( "gamma_1" = 2, "gamma_2" = 200, "delta_1" = 0.6, "delta_2" = 0.7, 
   "rho_1" = 0.5, "rho" = 0.6 )
bVrs <- c( b, "nu" = 1.1 )
xNames <- c( "X1", "X2", "X3" )


## checking cesCalc()
MishraCES$Y3 <- cesCalc( xNames = xNames, 
   data = MishraCES, coef = b, nested = TRUE )
MishraCES$Y3

# VRS
MishraCES$Y3Vrs <- cesCalc( xNames = xNames, 
   data = MishraCES, coef = bVrs, nested = TRUE )
MishraCES$Y3Vrs


