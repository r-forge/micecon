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
summary( cesLandLabor )
coef( cesLandLabor ) 
vcov( cesLandLabor ) 
coef( summary( cesLandLabor ) )

# variable returns to scale
cesLandLaborVrs <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   vrs = TRUE )
print.default( cesLandLaborVrs )
print( cesLandLaborVrs )
summary( cesLandLaborVrs )
coef( cesLandLaborVrs )
vcov( cesLandLaborVrs )
coef( summary( cesLandLaborVrs ) )

# using the CG optimization method
cesLandLaborCg <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   method = "CG" )
print.default( cesLandLaborCg )
print( cesLandLaborCg )
summary( cesLandLaborCg )
coef( cesLandLaborCg )
vcov( cesLandLaborCg )
coef( summary( cesLandLaborCg ) )

# using the SANN optimization method
set.seed( 123 )
cesLandLaborSann <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   method = "SANN" )
print.default( cesLandLaborSann )
print( cesLandLaborSann )
summary( cesLandLaborSann )
coef( cesLandLaborSann )
vcov( cesLandLaborSann )
coef( summary( cesLandLaborSann ) )

# using the BFGS optimization method
cesLandLaborBfgs <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   method = "BFGS" )
print.default( cesLandLaborBfgs )
print( cesLandLaborBfgs )
summary( cesLandLaborBfgs )
coef( cesLandLaborBfgs )
vcov( cesLandLaborBfgs )
coef( summary( cesLandLaborBfgs ) )

# using the L-BFGS-B optimization method with constrained alpha
cesLandLaborBfgsCon <- cesEst( "qOutput", c( "land", "qLabor" ),
   germanFarms, method = "L-BFGS-B", lower = c( -Inf, 0, -Inf ),
   upper = c( Inf, 1, Inf ) )
print.default( cesLandLaborBfgsCon )
print( cesLandLaborBfgsCon )
summary( cesLandLaborBfgsCon )
coef( cesLandLaborBfgsCon )
vcov( cesLandLaborBfgsCon )
coef( summary( cesLandLaborBfgsCon ) )

# Kmenta approximation
cesLandLaborKmenta <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
   vrs = TRUE, method = "Kmenta" )
print.default( cesLandLaborKmenta )
print( cesLandLaborKmenta )
summary( cesLandLaborKmenta )
coef( cesLandLaborKmenta )
vcov( cesLandLaborKmenta )
coef( summary( cesLandLaborKmenta ) )


## CES: Land & Intermediate Inputs
cesLandInt <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms )
print.default( cesLandInt )
print( cesLandInt )
summary( cesLandInt )
coef( summary( cesLandInt ) )

# variable returns to scale
cesLandIntVrs <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   vrs = TRUE )
print.default( cesLandIntVrs )
print( cesLandIntVrs )
summary( cesLandIntVrs )
coef( summary( cesLandIntVrs ) )

# using the CG optimization method
cesLandIntCg <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "CG" )
print.default( cesLandIntCg )
print( cesLandIntCg )
summary( cesLandIntCg )
coef( summary( cesLandIntCg ) )

# using the SANN optimization method
set.seed( 234 )
cesLandIntSann <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "SANN" )
print.default( cesLandIntSann )
print( cesLandIntSann )
summary( cesLandIntSann )
coef( summary( cesLandIntSann ) )

# using the BFGS optimization method
cesLandIntBfgs <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "BFGS" )
print.default( cesLandIntBfgs )
print( cesLandIntBfgs )
summary( cesLandIntBfgs )
coef( summary( cesLandIntBfgs ) )

# using the L-BFGS-B optimization method with constrained alpha
cesLandIntBfgsCon <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "L-BFGS-B", lower = c( -Inf, 0, -Inf ),
   upper = c( Inf, 1, Inf ) )
print.default( cesLandIntBfgsCon )
print( cesLandIntBfgsCon )
summary( cesLandIntBfgsCon )
coef( summary( cesLandIntBfgsCon ) )

# Kmenta approximation
cesLandIntKmenta <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
   method = "Kmenta" )
print.default( cesLandIntKmenta )
print( cesLandIntKmenta )
summary( cesLandIntKmenta )
coef( summary( cesLandIntKmenta ) )


############  cesCalc  ################
outLandLabor <- cesCalc( c( "land", "qLabor" ), germanFarms,
   coef( cesLandLabor ) )
print( outLandLabor )
all.equal( outLandLabor, cesCalc( c( "land", "qLabor" ), germanFarms,
   coef( cesLandLabor )[ c( 2, 3, 1 ) ] ) )
all.equal( outLandLabor, cesCalc( c( "land", "qLabor" ), germanFarms,
   unname( coef( cesLandLabor ) ) ) )

outLandLaborVrs <- cesCalc( c( "land", "qLabor" ), germanFarms,
   coef( cesLandLaborVrs ) )
print( outLandLaborVrs )
all.equal( outLandLaborVrs, cesCalc( c( "land", "qLabor" ), germanFarms,
   coef( cesLandLaborVrs )[ c( 3, 1, 4, 2 ) ] ) )
all.equal( outLandLaborVrs, cesCalc( c( "land", "qLabor" ), germanFarms,
   unname( coef( cesLandLaborVrs ) ) ) )
