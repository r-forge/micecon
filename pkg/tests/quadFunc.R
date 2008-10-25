library( micEcon )

## preparing data
data( germanFarms )
# output quantity:
germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
# quantity of variable inputs
germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
# a time trend to account for technical progress:
germanFarms$time <- c(1:20)

# estimate a quadratic production function
estResult <- quadFuncEst( "qOutput",
   c( "qLabor", "land", "qVarInput", "time" ), germanFarms )
coef( estResult )
print( estResult )

# compute fitted values
fitted <- quadFuncCalc( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ) )
all.equal( fitted, estResult$fitted )

# compute the marginal products of the inputs
margProducts <- quadFuncDeriv(
   c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ), vcov( estResult ) )
print( margProducts )
