library( micEcon )

# load data
data( germanFarms )
# output quantity:
germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
# quantity of intermediate inputs
germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput


## CES: Land & Labor
cesLandLabor <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms )
print( cesLandLabor )

# variable returns to scale
cesLandLaborVrs <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   vrs = TRUE )
print( cesLandLaborVrs )

# using the BFGS optimization method
cesLandLaborBfgs <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   method = "BFGS" )
print( cesLandLaborBfgs )


## CES: Land & Intermediate Inputs
cesLandInt <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms )
print( cesLandInt )

# variable returns to scale
cesLandIntVrs <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   vrs = TRUE )
print( cesLandIntVrs )

# using the BFGS optimization method
cesLandIntBfgs <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "BFGS" )
print( cesLandIntBfgs )
