library( micEconCES )

# load data
data( germanFarms )
# output quantity:
germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
# quantity of intermediate inputs
germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput


############  cesEstGridRho  ################
## CES: Land & Intermediate inputs
## Nelder-Mead, CRS
cesGridNm <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, rhoValues = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridNm ) 
print( cesGridNm )
summary( cesGridNm )
coef( cesGridNm )
vcov( cesGridNm )
coef( summary( cesGridNm ) )
fitted( cesGridNm )
residuals( cesGridNm )
plot( cesGridNm )

## Nelder-Mead, VRS
cesGridNmVrs <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, vrs = TRUE,
   rhoValues = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridNmVrs )
print( cesGridNmVrs )
summary( cesGridNmVrs )
coef( cesGridNmVrs )
vcov( cesGridNmVrs )
coef( summary( cesGridNmVrs ) )
plot( cesGridNmVrs )

# using the CG optimization method
cesGridCg <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "CG",
   rhoValues = seq( from = -0.8, to = 1.2, by = 0.4 ) )
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
cesGridCgVrs <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "CG", vrs = TRUE,
   rhoValues = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridCgVrs )
print( cesGridCgVrs )
summary( cesGridCgVrs )
coef( cesGridCgVrs )
vcov( cesGridCgVrs )
coef( summary( cesGridCgVrs ) )
plot( cesGridCgVrs )

# using the SANN optimization method
set.seed( 1 )
cesGridSann <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "SANN",
   rhoValues = seq( from = -0.8, to = 0.8, by = 0.8 ) )
print.default( cesGridSann )
print( cesGridSann )
summary( cesGridSann )
coef( cesGridSann )
vcov( cesGridSann )
coef( summary( cesGridSann ) )
plot( cesGridSann )

# using the SANN optimization method, VRS
set.seed( 21 )
cesGridSannVrs <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "SANN", vrs = TRUE,
   rhoValues = seq( from = -0.8, to = 0.8, by = 0.8 ) )
print.default( cesGridSannVrs )
print( cesGridSannVrs )
summary( cesGridSannVrs )
coef( cesGridSannVrs )
vcov( cesGridSannVrs )
coef( summary( cesGridSannVrs ) )
plot( cesGridSannVrs )

# using the BFGS optimization method
cesGridBfgs <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "BFGS",
   rhoValues = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridBfgs )
print( cesGridBfgs )
summary( cesGridBfgs )
coef( cesGridBfgs )
vcov( cesGridBfgs )
coef( summary( cesGridBfgs ) )
plot( cesGridBfgs )

# using the BFGS optimization method, VRS
cesGridBfgsVrs <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "BFGS", vrs = TRUE,
   rhoValues = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridBfgsVrs )
print( cesGridBfgsVrs )
summary( cesGridBfgsVrs )
coef( cesGridBfgsVrs )
vcov( cesGridBfgsVrs )
coef( summary( cesGridNmVrs ) )
plot( cesGridBfgsVrs )

# using the L-BFGS-B optimization method with constrained parameters
cesGridBfgsCon <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "L-BFGS-B",
   rhoValues = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridBfgsCon )
print( cesGridBfgsCon )
summary( cesGridBfgsCon )
coef( cesGridBfgsCon )
vcov( cesGridBfgsCon )
coef( summary( cesGridBfgsCon ) )
plot( cesGridBfgsCon )

# using the L-BFGS-B optimization method with constrained parameters, VRS
cesGridBfgsConVrs <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "L-BFGS-B", vrs = TRUE,
   rhoValues = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridBfgsConVrs )
print( cesGridBfgsConVrs )
summary( cesGridBfgsConVrs )
coef( cesGridBfgsConVrs )
vcov( cesGridBfgsConVrs )
coef( summary( cesGridBfgsConVrs ) )
plot( cesGridBfgsConVrs )

# using the Levenberg-Marquardt optimization method
cesGridLm <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "LM",
   rhoValues = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridLm )
print( cesGridLm )
summary( cesGridLm )
coef( cesGridLm )
vcov( cesGridLm )
coef( summary( cesGridLm ) )
plot( cesGridLm )

# using the Levenberg-Marquardt optimization method, VRS
cesGridLmVrs <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "LM", vrs = TRUE,
   rhoValues = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridLmVrs )
print( cesGridLmVrs )
summary( cesGridLmVrs )
coef( cesGridLmVrs )
vcov( cesGridLmVrs )
coef( summary( cesGridLmVrs ) )
plot( cesGridLmVrs )

# using the Newton-type optimization method implemented in nlm()
cesGridNewton <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "Newton",
   rhoValues = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridNewton )
print( cesGridNewton )
summary( cesGridNewton )
coef( cesGridNewton )
vcov( cesGridNewton )
coef( summary( cesGridNewton ) )
plot( cesGridNewton )

# using the Newton-type optimization method implemented in nlm(), VRS
cesGridNewtonVrs <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "Newton", vrs = TRUE,
   rhoValues = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridNewtonVrs )
print( cesGridNewtonVrs )
summary( cesGridNewtonVrs )
coef( cesGridNewtonVrs )
vcov( cesGridNewtonVrs )
coef( summary( cesGridNewtonVrs ) )
plot( cesGridNewtonVrs )

# using the PORT optimization rountine implemented in nlminb(), constrained
cesGridPort <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "PORT",
   rhoValues = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridPort )
print( cesGridPort )
summary( cesGridPort )
coef( cesGridPort )
vcov( cesGridNm )
coef( summary( cesGridPort ) )
plot( cesGridNm )

# using the PORT optimization rountine implemented in nlminb(), VRS, constrained
cesGridPortVrs <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "PORT", vrs = TRUE,
   rhoValues = seq( from = -0.8, to = 1.2, by = 0.4 ) )
print.default( cesGridPortVrs )
print( cesGridPortVrs )
summary( cesGridPortVrs )
coef( cesGridPortVrs )
vcov( cesGridPortVrs )
coef( summary( cesGridPortVrs ) )
plot( cesGridPortVrs )

# using the DE optimization method implemented in DEoptim(), CRS
set.seed( 321 )
cesGridDe <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "DE",
   rhoValues = seq( from = -0.8, to = 0.8, by = 0.8 ) )
print.default( cesGridDe )
print( cesGridDe )
summary( cesGridDe )
coef( cesGridDe )
vcov( cesGridDe )
coef( summary( cesGridDe ) )
plot( cesGridDe )

# using the DE optimization method implemented in DEoptim(), VRS
set.seed( 4321 )
cesGridDeVrs <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "DE", vrs = TRUE,
   rhoValues = seq( from = -0.8, to = 0.8, by = 0.8 ) )
print.default( cesGridDeVrs )
print( cesGridDeVrs )
summary( cesGridDeVrs )
coef( cesGridDeVrs )
vcov( cesGridDeVrs )
coef( summary( cesGridDeVrs ) )
plot( cesGridDeVrs )
