library( micEcon )

data( germanFarms )
# output quantity:
germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
# quantity of variable inputs
germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
# a time trend to account for technical progress:
germanFarms$time <- c(1:20)

# estimate a Cobb-Douglas production function
estResult <- translogEst( "qOutput", c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, linear = TRUE )

# calculate fitted values
fitted <- cobbDouglasCalc( c( "qLabor", "land", "qVarInput", "time" ),
   data = germanFarms, coef = coef( estResult )[ 1:5 ] )
print( fitted )
all.equal( fitted, estResult$fitted )

# calculate fitted values using logged independent variables
germanFarms$lQLabor    <- log( germanFarms$qLabor )
germanFarms$lLand      <- log( germanFarms$land )
germanFarms$lQVarInput <- log( germanFarms$qVarInput )
germanFarms$lTime      <- log( germanFarms$time )
fittedLogged <- cobbDouglasCalc( c( "lQLabor", "lLand", "lQVarInput", "lTime" ),
   data = germanFarms, coef = coef( estResult )[ 1:5 ], dataLogged = TRUE )
all.equal( fitted, exp( fittedLogged ) )

# coefficients not named
coefNoNames <- coef( estResult )[ 1:5 ]
names( coefNoNames ) <- NULL
fittedNoNames <- cobbDouglasCalc( c( "qLabor", "land", "qVarInput", "time" ),
   data = germanFarms, coef = coefNoNames )
all.equal( fitted, fittedNoNames )

# coefficients in a different order
coefDiffOrder <- coef( estResult )[ c( 3, 5, 1, 2, 4 ) ]
fittedDiffOrder <- cobbDouglasCalc( c( "qLabor", "land", "qVarInput", "time" ),
   data = germanFarms, coef = coefDiffOrder )
all.equal( fitted, fittedDiffOrder )
