library( micEcon )

## preparing data
data( germanFarms )
# output quantity:
germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
# quantity of variable inputs
germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
# a time trend to account for technical progress:
germanFarms$time <- c(1:20)

## testing translogEst
# estimate a translog production function
estResult <- translogEst( "qOutput", c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms )

print( estResult )

## testing translogCalc
fitted <- translogCalc( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$allCoef )

all.equal( fitted, estResult$fitted )

## testing translogDeriv
margProducts <- translogDeriv( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$allCoef, estResult$allCoefCov )

print( margProducts )

## testing translogHessian
# compute the Hessian matrices
hessians <- translogHessian( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$allCoef )

print( hessians )

# compute the bordered Hessian matrices
borderedHessians <- translogHessian( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$allCoef, bordered = TRUE )

print( borderedHessians )

