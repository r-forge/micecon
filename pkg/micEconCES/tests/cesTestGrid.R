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
   data = germanFarms, by = 0.4, to = 1.2 )
print.default( cesGridNm ) 
print( cesGridNm )
summary( cesGridNm )
coef( cesGridNm )
vcov( cesGridNm )
coef( summary( cesGridNm ) )
fitted( cesGridNm )
residuals( cesGridNm )

## Nelder-Mead, VRS
cesGridNmVrs <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, vrs = TRUE, by = 0.4, to = 1.2 )
print.default( cesGridNmVrs )
print( cesGridNmVrs )
summary( cesGridNmVrs )
coef( cesGridNmVrs )
vcov( cesGridNmVrs )
coef( summary( cesGridNmVrs ) )

# using the CG optimization method
cesGridCg <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "CG", by = 0.4, to = 1.2 )
print.default( cesGridCg )
print( cesGridCg )
summary( cesGridCg )
coef( cesGridCg )
vcov( cesGridCg )
coef( summary( cesGridCg ) )
fitted( cesGridCg )
residuals( cesGridCg )

# using the CG optimization method, VRS
cesGridCgVrs <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "CG", vrs = TRUE, by = 0.4, to = 1.2 )
print.default( cesGridCgVrs )
print( cesGridCgVrs )
summary( cesGridCgVrs )
coef( cesGridCgVrs )
vcov( cesGridCgVrs )
coef( summary( cesGridCgVrs ) )

# using the SANN optimization method
set.seed( 1 )
cesGridSann <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "SANN", by = 0.8, to = 0.8 )
print.default( cesGridSann )
print( cesGridSann )
summary( cesGridSann )
coef( cesGridSann )
vcov( cesGridSann )
coef( summary( cesGridSann ) )

# using the SANN optimization method, VRS
set.seed( 21 )
cesGridSannVrs <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "SANN", vrs = TRUE, by = 0.8, to = 0.8 )
print.default( cesGridSannVrs )
print( cesGridSannVrs )
summary( cesGridSannVrs )
coef( cesGridSannVrs )
vcov( cesGridSannVrs )
coef( summary( cesGridSannVrs ) )

# using the BFGS optimization method
cesGridBfgs <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "BFGS", by = 0.4, to = 1.2 )
print.default( cesGridBfgs )
print( cesGridBfgs )
summary( cesGridBfgs )
coef( cesGridBfgs )
vcov( cesGridBfgs )
coef( summary( cesGridBfgs ) )

# using the BFGS optimization method, VRS
cesGridBfgsVrs <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "BFGS", vrs = TRUE, by = 0.4, to = 1.2 )
print.default( cesGridBfgsVrs )
print( cesGridBfgsVrs )
summary( cesGridBfgsVrs )
coef( cesGridBfgsVrs )
vcov( cesGridBfgsVrs )
coef( summary( cesGridNmVrs ) )

# using the L-BFGS-B optimization method with constrained parameters
cesGridBfgsCon <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "L-BFGS-B", by = 0.4, to = 1.2 )
print.default( cesGridBfgsCon )
print( cesGridBfgsCon )
summary( cesGridBfgsCon )
coef( cesGridBfgsCon )
vcov( cesGridBfgsCon )
coef( summary( cesGridBfgsCon ) )

# using the L-BFGS-B optimization method with constrained parameters, VRS
cesGridBfgsConVrs <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "L-BFGS-B", vrs = TRUE, by = 0.4, to = 1.2 )
print.default( cesGridBfgsConVrs )
print( cesGridBfgsConVrs )
summary( cesGridBfgsConVrs )
coef( cesGridBfgsConVrs )
vcov( cesGridBfgsConVrs )
coef( summary( cesGridBfgsConVrs ) )

# using the Levenberg-Marquardt optimization method
cesGridLm <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "LM", by = 0.4, to = 1.2 )
print.default( cesGridLm )
print( cesGridLm )
summary( cesGridLm )
coef( cesGridLm )
vcov( cesGridLm )
coef( summary( cesGridLm ) )

# using the Levenberg-Marquardt optimization method, VRS
cesGridLmVrs <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "LM", vrs = TRUE, by = 0.4, to = 1.2 )
print.default( cesGridLmVrs )
print( cesGridLmVrs )
summary( cesGridLmVrs )
coef( cesGridLmVrs )
vcov( cesGridLmVrs )
coef( summary( cesGridLmVrs ) )

# using the Newton-type optimization method implemented in nlm()
cesGridNewton <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "Newton", by = 0.4, to = 1.2 )
print.default( cesGridNewton )
print( cesGridNewton )
summary( cesGridNewton )
coef( cesGridNewton )
vcov( cesGridNewton )
coef( summary( cesGridNewton ) )

# using the Newton-type optimization method implemented in nlm(), VRS
cesGridNewtonVrs <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "Newton", vrs = TRUE, by = 0.4, to = 1.2 )
print.default( cesGridNewtonVrs )
print( cesGridNewtonVrs )
summary( cesGridNewtonVrs )
coef( cesGridNewtonVrs )
vcov( cesGridNewtonVrs )
coef( summary( cesGridNewtonVrs ) )

# using the PORT optimization rountine implemented in nlminb(), constrained
cesGridPort <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "PORT", by = 0.4, to = 1.2 )
print.default( cesGridPort )
print( cesGridPort )
summary( cesGridPort )
coef( cesGridPort )
vcov( cesGridNm )
coef( summary( cesGridPort ) )

# using the PORT optimization rountine implemented in nlminb(), VRS, constrained
cesGridPortVrs <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "PORT", vrs = TRUE, by = 0.4, to = 1.2 )
print.default( cesGridPortVrs )
print( cesGridPortVrs )
summary( cesGridPortVrs )
coef( cesGridPortVrs )
vcov( cesGridPortVrs )
coef( summary( cesGridPortVrs ) )

# using the DE optimization method implemented in DEoptim(), CRS
set.seed( 321 )
cesGridDe <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "DE", by = 0.8, to = 0.8 )
print.default( cesGridDe )
print( cesGridDe )
summary( cesGridDe )
coef( cesGridDe )
vcov( cesGridDe )
coef( summary( cesGridDe ) )

# using the DE optimization method implemented in DEoptim(), VRS
set.seed( 4321 )
cesGridDeVrs <- cesEstGridRho( yName = "qOutput", xNames = c( "land", "qVarInput" ),
   data = germanFarms, method = "DE", vrs = TRUE, by = 0.8, to = 0.8 )
print.default( cesGridDeVrs )
print( cesGridDeVrs )
summary( cesGridDeVrs )
coef( cesGridDeVrs )
vcov( cesGridDeVrs )
coef( summary( cesGridDeVrs ) )
