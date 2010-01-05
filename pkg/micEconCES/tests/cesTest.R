library( micEconCES )

# load data
data( germanFarms )
# output quantity:
germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
# quantity of intermediate inputs
germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput


## CES: Land & Labor
cesLandLabor <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms )
print.default( cesLandLabor ) 
print( cesLandLabor )
coef( cesLandLabor ) 
vcov( cesLandLabor ) 
summary( cesLandLabor )$coefficients

# variable returns to scale
cesLandLaborVrs <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   vrs = TRUE )
print.default( cesLandLaborVrs )
print( cesLandLaborVrs )
coef( cesLandLaborVrs )
vcov( cesLandLaborVrs )
summary( cesLandLaborVrs )$coefficients

# using the BFGS optimization method
cesLandLaborBfgs <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   method = "BFGS" )
print.default( cesLandLaborBfgs )
print( cesLandLaborBfgs )
coef( cesLandLaborBfgs )
vcov( cesLandLaborBfgs )
summary( cesLandLaborBfgs )$coefficients

# using the L-BFGS-B optimization method with constrained alpha
cesLandLaborBfgsCon <- cesEst( "qOutput", c( "land", "qLabor" ),
   germanFarms, method = "L-BFGS-B", lower = c( -Inf, 0, -Inf ),
   upper = c( Inf, 1, Inf ) )
print.default( cesLandLaborBfgsCon )
print( cesLandLaborBfgsCon )
coef( cesLandLaborBfgsCon )
vcov( cesLandLaborBfgsCon )
summary( cesLandLaborBfgsCon )$coefficients


## CES: Land & Intermediate Inputs
cesLandInt <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms )
print.default( cesLandInt )
print( cesLandInt )
summary( cesLandInt )$coefficients

# variable returns to scale
cesLandIntVrs <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   vrs = TRUE )
print.default( cesLandIntVrs )
print( cesLandIntVrs )
summary( cesLandIntVrs )$coefficients

# using the BFGS optimization method
cesLandIntBfgs <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "BFGS" )
print.default( cesLandIntBfgs )
print( cesLandIntBfgs )
summary( cesLandIntBfgs )$coefficients

# using the L-BFGS-B optimization method with constrained alpha
cesLandIntBfgsCon <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "L-BFGS-B", lower = c( -Inf, 0, -Inf ),
   upper = c( Inf, 1, Inf ) )
print.default( cesLandIntBfgsCon )
print( cesLandIntBfgsCon )
summary( cesLandIntBfgsCon )$coefficients
