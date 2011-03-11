library( "micEconCES" )
data( "GermanIndustry" )

print( GermanIndustry )

GermanIndustry$time <- GermanIndustry$year - 1963

b <- c( "gamma" = 38, "lambda" = 0.0222, "delta_1" = 0.5, "delta_2" = 0.1, 
   "rho_1" = 0.5300, "rho" = 0.1813 )

GermanIndustry$YCalc <- cesCalc( xNames = c( "K", "E", "A" ), tName = "time",
   data = GermanIndustry, coef = b, nested = TRUE )

GermanIndustry$YCalc 


