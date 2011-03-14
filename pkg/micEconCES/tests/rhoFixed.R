# load the micEconCES package
library( micEconCES )

# load data
data( germanFarms )
# output quantity:
germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
# quantity of intermediate inputs
germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
# time trend
germanFarms$time <- c( 1:nrow( germanFarms ) )

# names of explanatory variables
xxNames <- c( "land", "qVarInput" )


## Nelder-Mead, CRS
cesNm <- cesEst( "qOutput", xxNames, germanFarms, rho = 2, method = "NM" )
print.default( cesNm ) 
print( cesNm )
summary( cesNm )
coef( cesNm ) 
vcov( cesNm ) 
coef( summary( cesNm ) )

## Nelder-Mead, VRS
cesNmVrs <- cesEst( "qOutput", xxNames, germanFarms, vrs = TRUE, rho = -0.1,
   method = "Nelder-Mead" )
print.default( cesNmVrs )
print( cesNmVrs )
summary( cesNmVrs )
coef( cesNmVrs )
vcov( cesNmVrs )
coef( summary( cesNmVrs ) )

## Nelder-Mead, TC, VRS
cesNmTcVrs <- cesEst( "qOutput", xxNames, tName = "time", data = germanFarms, 
   vrs = TRUE, rho = -0.4, method = "Nelder-Mead" )
print.default( cesNmTcVrs )
print( cesNmTcVrs )
summary( cesNmTcVrs )

## Conjugate Gradients, CRS
cesCg <- cesEst( "qOutput", xxNames, germanFarms, method = "CG", rho = 2 )
print.default( cesCg )
print( cesCg )
summary( cesCg )
coef( cesCg )
vcov( cesCg )
coef( summary( cesCg ) )

## Conjugate Gradients, TC, CRS
cesCgTc <- cesEst( "qOutput", xxNames, tName = "time", data = germanFarms, 
   method = "CG", rho = 0.2 )
print.default( cesCgTc )
print( cesCgTc )
summary( cesCgTc )

## Conjugate Gradients, VRS
cesCgVrs <- cesEst( "qOutput", xxNames, germanFarms, method = "CG", vrs = TRUE,
   rho = -0.1 )
print.default( cesCgVrs )
print( cesCgVrs )
summary( cesCgVrs )
coef( cesCgVrs )
vcov( cesCgVrs )
coef( summary( cesCgVrs ) )

## Simulated Annealing, CRS
cesSann <- cesEst( "qOutput", xxNames, germanFarms, method = "SANN", rho = 2 )
print.default( cesSann )
print( cesSann )
summary( cesSann )
coef( cesSann )
vcov( cesSann )
coef( summary( cesSann ) )

## Simulated Annealing, TC, CRS
cesSannTc <- cesEst( "qOutput", xxNames, tName = "time", data = germanFarms, 
   method = "SANN", rho = 0 )
print.default( cesSannTc )
print( cesSannTc )
summary( cesSannTc )

## Simulated Annealing, VRS
cesSannVrs <- cesEst( "qOutput", xxNames, germanFarms, method = "SANN", vrs = TRUE,
   rho = -0.1 )
print.default( cesSannVrs )
print( cesSannVrs )
summary( cesSannVrs )
coef( cesSannVrs )
vcov( cesSannVrs )
coef( summary( cesSannVrs ) )

## BFGS, CRS
cesBfgs <- cesEst( "qOutput", xxNames, germanFarms, method = "BFGS", rho = 2 )
print.default( cesBfgs )
print( cesBfgs )
summary( cesBfgs )
coef( cesBfgs )
vcov( cesBfgs )
coef( summary( cesBfgs ) )

## BFGS, VRS
cesBfgsVrs <- cesEst( "qOutput", xxNames, germanFarms, method = "BFGS", vrs = TRUE,
   rho = -0.1 )
print.default( cesBfgsVrs )
print( cesBfgsVrs )
summary( cesBfgsVrs )
coef( cesBfgsVrs )
vcov( cesBfgsVrs )
coef( summary( cesBfgsVrs ) )

## BFGS, TC, VRS
cesBfgsTcVrs <- cesEst( "qOutput", xxNames, tName = "time", 
   data = germanFarms, method = "BFGS", vrs = TRUE, rho = -0.4 )
print.default( cesBfgsTcVrs )
print( cesBfgsTcVrs )
summary( cesBfgsTcVrs )

## L-BFGS-B with constrained parameters, CRS
cesBfgsCon <- cesEst( "qOutput", xxNames, germanFarms, method = "L-BFGS-B",
   rho = 2 )
print.default( cesBfgsCon )
print( cesBfgsCon )
summary( cesBfgsCon )
coef( cesBfgsCon )
vcov( cesBfgsCon )
coef( summary( cesBfgsCon ) )

## L-BFGS-B with constrained parameters, TC, CRS
cesBfgsConTc <- cesEst( "qOutput", xxNames, tName = "time", data = germanFarms, 
   method = "L-BFGS-B", rho = 1.2 )
print.default( cesBfgsConTc )
print( cesBfgsConTc )
summary( cesBfgsConTc )

## L-BFGS-B with constrained parameters, VRS
cesBfgsConVrs <- cesEst( "qOutput", xxNames, germanFarms, method = "L-BFGS-B",
   vrs = TRUE, rho = -0.1 )
print.default( cesBfgsConVrs )
print( cesBfgsConVrs )
summary( cesBfgsConVrs )
coef( cesBfgsConVrs )
vcov( cesBfgsConVrs )
coef( summary( cesBfgsConVrs ) )

## Levenberg-Marquardt, CRS
cesLm <- cesEst( "qOutput", xxNames, germanFarms,
   control = nls.lm.control( maxiter = 200 ), rho = 2 )
print.default( cesLm )
print( cesLm )
summary( cesLm )
coef( cesLm )
vcov( cesLm )
coef( summary( cesLm ) )

## Levenberg-Marquardt, VRS
cesLmVrs <- cesEst( "qOutput", xxNames, germanFarms, vrs = TRUE,
   control = nls.lm.control( maxiter = 200 ), rho = -0.1 )
print.default( cesLmVrs )
print( cesLmVrs )
summary( cesLmVrs )
coef( cesLmVrs )
vcov( cesLmVrs )
coef( summary( cesLmVrs ) )

## Levenberg-Marquardt, TC, VRS
cesLmTcVrs <- cesEst( "qOutput", xxNames, tName = "time", germanFarms, 
   vrs = TRUE, control = nls.lm.control( maxiter = 200 ), rho = 0.7 )
print.default( cesLmTcVrs )
print( cesLmTcVrs )
summary( cesLmTcVrs )

## Newton-type, CRS
cesNewton <- cesEst( "qOutput", xxNames, germanFarms, method = "Newton", rho = 2 )
print.default( cesNewton )
print( cesNewton )
summary( cesNewton )
coef( cesNewton )
vcov( cesNewton )
coef( summary( cesNewton ) )

## Newton-type, TC, CRS
cesNewtonTc <- cesEst( "qOutput", xxNames, tName = "time", data = germanFarms, 
   method = "Newton", rho = 0.8 )
print.default( cesNewtonTc )
print( cesNewtonTc )
summary( cesNewtonTc )

## Newton-type, VRS
cesNewtonVrs <- cesEst( "qOutput", xxNames, germanFarms, method = "Newton",
   vrs = TRUE, rho = -0.1 )
print.default( cesNewtonVrs )
print( cesNewtonVrs )
summary( cesNewtonVrs )
coef( cesNewtonVrs )
vcov( cesNewtonVrs )
coef( summary( cesNewtonVrs ) )

## PORT, CRS
cesPort <- cesEst( "qOutput", xxNames, germanFarms, method = "PORT", rho = 2 )
print.default( cesPort )
print( cesPort )
summary( cesPort )
coef( cesPort )
vcov( cesPort )
coef( summary( cesPort ) )

## PORT, VRS
cesPortVrs <- cesEst( "qOutput", xxNames, germanFarms, method = "PORT", vrs = TRUE,
   rho = -0.1 )
print.default( cesPortVrs )
print( cesPortVrs )
summary( cesPortVrs )
coef( cesPortVrs )
vcov( cesPortVrs )
coef( summary( cesPortVrs ) )

## PORT, TC, VRS
cesPortTcVrs <- cesEst( "qOutput", xxNames, tName = "time", data = germanFarms, 
   method = "PORT", vrs = TRUE, rho = 0 )
print.default( cesPortTcVrs )
print( cesPortTcVrs )
summary( cesPortTcVrs )

## DE, CRS
cesDe <- cesEst( "qOutput", xxNames, germanFarms, method = "DE",
   control = DEoptim.control( trace = FALSE ), rho = 2 )
print.default( cesDe )
print( cesDe )
summary( cesDe )
coef( cesDe )
vcov( cesDe )
coef( summary( cesDe ) )

## DE, TC, CRS
cesDeTc <- cesEst( "qOutput", xxNames, tName = "time", data = germanFarms, 
   method = "DE", control = DEoptim.control( trace = FALSE ), rho = 0.5 )
print.default( cesDeTc )
print( cesDeTc )
summary( cesDeTc )

## DE, VRS
cesDeVrs <- cesEst( "qOutput", xxNames, germanFarms, method = "DE", vrs = TRUE,
   control = DEoptim.control( trace = FALSE ), rho = -0.1, random.seed = 234 )
print.default( cesDeVrs )
print( cesDeVrs )
summary( cesDeVrs )
coef( cesDeVrs )
vcov( cesDeVrs )
coef( summary( cesDeVrs ) )

