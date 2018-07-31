library( "micEconCES" )
options( digits = 3 )
set.seed( 12345 )  # to make the bootstrapping in dwt() reproducible

# load data
data( germanFarms, package = "micEcon" )
# output quantity:
germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
# quantity of intermediate inputs
germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
# time trend
germanFarms$time <- c( 1:nrow( germanFarms ) )


############  cesEstGridRho  ################
## CES: Land & Intermediate inputs
## Nelder-Mead, CRS
cesGridNm <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, rho = seq( from = -0.8, to = 1.2, by = 0.4 ),
   method = "Nelder-Mead", returnGrad = TRUE )
print.default( cesGridNm ) 
print( cesGridNm )
summary( cesGridNm )
coef( cesGridNm )
vcov( cesGridNm )
coef( summary( cesGridNm ) )
fitted( cesGridNm )
residuals( cesGridNm )
plot( cesGridNm )
dwt( cesGridNm )

## Nelder-Mead, TC, CRS
cesGridNmTc <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   tName = "time", data = germanFarms, method = "Nelder-Mead", 
   rho = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridNmTc ) 
print( cesGridNmTc )
summary( cesGridNmTc )
try( dwt( summary( cesGridNmTc ) ) )

## Nelder-Mead, VRS
cesGridNmVrs <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "NM", vrs = TRUE,
   rho = seq( from = -0.8, to = 1.2, by = 0.4 ), returnGrad = TRUE )
print.default( cesGridNmVrs )
print( cesGridNmVrs )
summary( cesGridNmVrs )
coef( cesGridNmVrs )
vcov( cesGridNmVrs )
coef( summary( cesGridNmVrs ) )
plot( cesGridNmVrs )
dwt( cesGridNmVrs )

# using the CG optimization method
cesGridCg <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "CG",
   rho = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridCg )
print( cesGridCg )
summary( cesGridCg )
coef( cesGridCg )
vcov( cesGridCg )
coef( summary( cesGridCg ) )
fitted( cesGridCg )
residuals( cesGridCg )
plot( cesGridCg )

# using the CG optimization method, VRS
cesGridCgVrs <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "CG", vrs = TRUE,
   rho = seq( from = -0.8, to = 1.2, by = 0.4 ), returnGrad = TRUE )
print.default( cesGridCgVrs )
print( cesGridCgVrs )
summary( cesGridCgVrs )
coef( cesGridCgVrs )
vcov( cesGridCgVrs )
coef( summary( cesGridCgVrs ) )
plot( cesGridCgVrs )
dwt( cesGridCgVrs )

# using the CG optimization method, TC, VRS
cesGridCgTcVrs <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   tName = "time", data = germanFarms, method = "CG", vrs = TRUE,
   rho = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridCgTcVrs )
print( cesGridCgTcVrs )
summary( cesGridCgTcVrs )

# using the SANN optimization method
cesGridSann <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "SANN", random.seed = 1,
   rho = seq( from = -0.8, to = 0.8, by = 0.8 ) )
print.default( cesGridSann )
print( cesGridSann )
summary( cesGridSann )
coef( cesGridSann )
vcov( cesGridSann )
coef( summary( cesGridSann ) )
plot( cesGridSann )

# using the SANN optimization method, TC, CRS
cesGridSannTc <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   tName = "time", data = germanFarms, method = "SANN", random.seed = 31,
   rho = seq( from = -0.8, to = 0.8, by = 0.8 ) )
print.default( cesGridSannTc )
print( cesGridSannTc )
summary( cesGridSannTc )

# using the SANN optimization method, VRS
cesGridSannVrs <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "SANN", vrs = TRUE, random.seed = 21,
   rho = seq( from = -0.8, to = 0.8, by = 0.8 ) )
print.default( cesGridSannVrs )
print( cesGridSannVrs )
summary( cesGridSannVrs )
coef( cesGridSannVrs )
vcov( cesGridSannVrs )
coef( summary( cesGridSannVrs ) )
plot( cesGridSannVrs )

# using the BFGS optimization method
cesGridBfgs <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "BFGS",
   rho = seq( from = -0.8, to = 1.2, by = 0.4 ), returnGrad = TRUE )
print.default( cesGridBfgs )
print( cesGridBfgs )
summary( cesGridBfgs )
coef( cesGridBfgs )
vcov( cesGridBfgs )
coef( summary( cesGridBfgs ) )
plot( cesGridBfgs )
dwt( cesGridBfgs )

# using the BFGS optimization method, VRS
cesGridBfgsVrs <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "BFGS", vrs = TRUE,
   rho = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridBfgsVrs )
print( cesGridBfgsVrs )
summary( cesGridBfgsVrs )
coef( cesGridBfgsVrs )
vcov( cesGridBfgsVrs )
coef( summary( cesGridNmVrs ) )
plot( cesGridBfgsVrs )

# using the BFGS optimization method, TC, VRS
cesGridBfgsTcVrs <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   tName = "time", data = germanFarms, method = "BFGS", vrs = TRUE,
   rho = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridBfgsTcVrs )
print( cesGridBfgsTcVrs )
summary( cesGridBfgsTcVrs )

# using the BFGS optimization method, TC, multErr, VRS
cesGridBfgsTcMeVrs <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   tName = "time", data = germanFarms, method = "BFGS", vrs = TRUE,
   multErr = TRUE, rho = seq( from = -0.4, to = 1.2, by = 0.4 ),
   control = list( maxit = 1000 ) )
print.default( cesGridBfgsTcMeVrs )
print( cesGridBfgsTcMeVrs )
summary( cesGridBfgsTcMeVrs )

# using the L-BFGS-B optimization method with constrained parameters
cesGridBfgsCon <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "L-BFGS-B",
   rho = seq( from = -0.8, to = 1.2, by = 0.4 ), returnGrad = TRUE )
print.default( cesGridBfgsCon )
print( cesGridBfgsCon )
summary( cesGridBfgsCon )
coef( cesGridBfgsCon )
vcov( cesGridBfgsCon )
coef( summary( cesGridBfgsCon ) )
plot( cesGridBfgsCon )
dwt( cesGridBfgsCon )

# using the L-BFGS-B optimization method, TC, CRS
cesGridBfgsConTc <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   tName = "time", data = germanFarms, method = "L-BFGS-B",
   rho = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridBfgsConTc )
print( cesGridBfgsConTc )
summary( cesGridBfgsConTc )

# using the L-BFGS-B optimization method with constrained parameters, VRS
cesGridBfgsConVrs <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "L-BFGS-B", vrs = TRUE,
   rho = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridBfgsConVrs )
print( cesGridBfgsConVrs )
summary( cesGridBfgsConVrs )
coef( cesGridBfgsConVrs )
vcov( cesGridBfgsConVrs )
coef( summary( cesGridBfgsConVrs ) )
plot( cesGridBfgsConVrs )

# using the Levenberg-Marquardt optimization method
cesGridLm <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, rho = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridLm )
print( cesGridLm )
summary( cesGridLm )
coef( cesGridLm )
vcov( cesGridLm )
coef( summary( cesGridLm ) )
plot( cesGridLm )

# using the Levenberg-Marquardt optimization method, VRS
cesGridLmVrs <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, vrs = TRUE,
   rho = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridLmVrs )
print( cesGridLmVrs )
summary( cesGridLmVrs )
coef( cesGridLmVrs )
vcov( cesGridLmVrs )
coef( summary( cesGridLmVrs ) )
plot( cesGridLmVrs )

# using the Levenberg-Marquardt optimization method, TC, VRS
cesGridLmTcVrs <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   tName = "time", data = germanFarms, vrs = TRUE,
   rho = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridLmTcVrs )
print( cesGridLmTcVrs )
summary( cesGridLmTcVrs )

# using the Levenberg-Marquardt optimization method, multErr, VRS
cesGridLmMeVrs <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, vrs = TRUE, multErr = TRUE,
   rho = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridLmMeVrs )
print( cesGridLmMeVrs )
summary( cesGridLmMeVrs )
vcov( cesGridLmMeVrs )

# using the Newton-type optimization method implemented in nlm()
cesGridNewton <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "Newton",
   rho = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridNewton )
print( cesGridNewton )
summary( cesGridNewton )
coef( cesGridNewton )
vcov( cesGridNewton )
coef( summary( cesGridNewton ) )
plot( cesGridNewton )

# using the Newton-type optimization method implemented in nlm(), TC, CRS
cesGridNewtonTc <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   tName = "time", data = germanFarms, method = "Newton",
   rho = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridNewtonTc )
print( cesGridNewtonTc )
summary( cesGridNewtonTc )

# using the Newton-type optimization method implemented in nlm(), VRS
cesGridNewtonVrs <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "Newton", vrs = TRUE,
   rho = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridNewtonVrs )
print( cesGridNewtonVrs )
summary( cesGridNewtonVrs )
coef( cesGridNewtonVrs )
vcov( cesGridNewtonVrs )
coef( summary( cesGridNewtonVrs ) )
plot( cesGridNewtonVrs )

# using the Newton-type optimization method implemented in nlm(), TC, multErr, CRS
cesGridNewtonTcMe <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   tName = "time", data = germanFarms, method = "Newton", multErr = TRUE,
   rho = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridNewtonTcMe )
print( cesGridNewtonTcMe )
summary( cesGridNewtonTcMe )

# using the PORT optimization rountine implemented in nlminb(), constrained
cesGridPort <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "PORT",
   rho = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridPort )
print( cesGridPort )
summary( cesGridPort )
coef( cesGridPort )
vcov( cesGridPort )
coef( summary( cesGridPort ) )
plot( cesGridPort )

# using the PORT optimization rountine implemented in nlminb(), VRS, constrained
cesGridPortVrs <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "PORT", vrs = TRUE,
   rho = seq( from = -0.8, to = 1.2, by = 0.4 ), returnGrad = TRUE )
print.default( cesGridPortVrs )
print( cesGridPortVrs )
summary( cesGridPortVrs )
coef( cesGridPortVrs )
vcov( cesGridPortVrs )
coef( summary( cesGridPortVrs ) )
plot( cesGridPortVrs )
dwt( cesGridPortVrs )

# using the PORT optimization rountine implemented in nlminb(), TC, VRS
cesGridPortTcVrs <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   tName = "time", data = germanFarms, method = "PORT", vrs = TRUE,
   rho = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridPortTcVrs )
print( cesGridPortTcVrs )
summary( cesGridPortTcVrs )

# using the PORT optimization rountine implemented in nlminb(), multErr
cesGridPortMe <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "PORT", multErr = TRUE,
   rho = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridPortMe )
print( cesGridPortMe )
summary( cesGridPortMe )
vcov( cesGridPortMe )

# using the DE optimization method implemented in DEoptim(), CRS
cesGridDe <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "DE", random.seed = 321,
   rho = seq( from = -0.8, to = 0.8, by = 0.8 ),
   control = DEoptim.control( trace = FALSE ) )
print.default( cesGridDe )
print( cesGridDe )
summary( cesGridDe )
coef( cesGridDe )
vcov( cesGridDe )
coef( summary( cesGridDe ) )
plot( cesGridDe )

# using the DE optimization method implemented in DEoptim(), TC, CRS
cesGridDeTc <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   tName = "time", data = germanFarms, method = "DE", random.seed = 321,
   rho = seq( from = -0.8, to = 0.8, by = 0.8 ), returnGrad = TRUE,
   control = DEoptim.control( trace = FALSE ) )
print.default( cesGridDeTc )
print( cesGridDeTc )
summary( cesGridDeTc )
dwt( cesGridDeTc )

# using the DE optimization method implemented in DEoptim(), VRS
cesGridDeVrs <- cesEst( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "DE", vrs = TRUE, random.seed = 4321,
   rho = seq( from = -0.8, to = 0.8, by = 0.8 ),
   control = DEoptim.control( trace = FALSE ) )
print.default( cesGridDeVrs )
print( cesGridDeVrs )
summary( cesGridDeVrs )
coef( cesGridDeVrs )
vcov( cesGridDeVrs )
coef( summary( cesGridDeVrs ) )
plot( cesGridDeVrs )
