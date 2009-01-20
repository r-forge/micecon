library( micEcon )
library( plm )

## preparing data
data( germanFarms )
# output quantity:
germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
# quantity of variable inputs
germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
# a time trend to account for technical progress:
germanFarms$time <- c(1:20)

## estimate a quadratic production function
estResult <- quadFuncEst( "qOutput",
   c( "qLabor", "land", "qVarInput", "time" ), germanFarms )
coef( estResult )
print( estResult )

## compute fitted values
fitted <- quadFuncCalc( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ) )
all.equal( fitted, estResult$fitted )

## compute the marginal products of the inputs
margProducts <- quadFuncDeriv(
   c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ), vcov( estResult ) )
print( margProducts )

## estimate a quadratic production function with a shifter
estResultShifter <- quadFuncEst( yName = "qOutput",
   xNames = c( "qLabor", "land", "qVarInput" ),
   shifterNames = "time", data = germanFarms )
coef( estResultShifter )
print( estResultShifter )
# compute fitted values
fitted <- quadFuncCalc( c( "qLabor", "land", "qVarInput" ),
   shifterNames = "time", data = germanFarms, coef( estResultShifter ) )
all.equal( fitted, estResultShifter$fitted )
# compute marginal products = partial derivatives
margProdShifter <- quadFuncDeriv(
   c( "qLabor", "land", "qVarInput" ),
   germanFarms, coef( estResultShifter ), vcov( estResultShifter ) )
print( margProdShifter )

## estimate a quadratic production function with 2 shifters
germanFarms$timeSq <- germanFarms$time^2
estResultShifter2 <- quadFuncEst( yName = "qOutput",
   xNames = c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "time", "timeSq" ), data = germanFarms )
coef( estResultShifter2 )
print( estResultShifter2 )
# compute fitted values
fitted <- quadFuncCalc( c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "time", "timeSq" ), data = germanFarms, 
   coef( estResultShifter2 ) )
all.equal( fitted, estResultShifter2$fitted )

## estimate a linear functions with quadFuncEst
estResultLinear <- quadFuncEst( yName = "qOutput", xNames = NULL,
   shifterNames = c( "time", "qLabor", "land", "qVarInput" ),
   data = germanFarms )
coef( estResultLinear )
print( estResultLinear )
estResultLin <- quadFuncEst( yName = "qOutput", 
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   data = germanFarms, linear = TRUE )
all.equal( coef( estResultLinear ), coef( estResultLin )[1:5],
   check.attributes = FALSE )
all.equal( vcov( estResultLinear ), vcov( estResultLin )[1:5,1:5],
   check.attributes = FALSE )
coef( estResultLin )
vcov( estResultLin )
print( estResultLin )
# compute fitted values
fitted <- quadFuncCalc( xNames = NULL,
   shifterNames = c( "time", "qLabor", "land", "qVarInput" ), 
   data = germanFarms, coef( estResultLinear ) )
all.equal( fitted, estResultLinear$fitted )
fitted <- quadFuncCalc( xNames = c( "time", "qLabor", "land", "qVarInput" ), 
   data = germanFarms, coef( estResultLin ) )
all.equal( fitted, estResultLin$fitted )
all.equal( estResultLinear$fitted, estResultLin$fitted )
# compute partial derivatives
margProducts <- quadFuncDeriv(
   c( "qLabor", "land", "qVarInput", "time" ),
   data = germanFarms, coef = coef( estResultLin ), 
   coefCov = vcov( estResultLin ) )
sd( margProducts$deriv )
all.equal( margProducts$deriv[1,], coef( estResultLin )[2:5], 
   check.attributes = FALSE )
all.equal( margProducts$variance[1,], diag( vcov( estResultLin ) )[2:5], 
   check.attributes = FALSE )

## estimate a quadratic production function with a logical variable as shifter
germanFarms$reUnif <- germanFarms$time >= 16
estResultShifterLogi <- quadFuncEst( yName = "qOutput",
   xNames = c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "reUnif" ), data = germanFarms )
coef( estResultShifterLogi )
print( estResultShifterLogi )
# compute fitted values
fitted <- quadFuncCalc( c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "reUnif" ), data = germanFarms, 
   coef( estResultShifterLogi ) )
all.equal( fitted, estResultShifterLogi$fitted )

## estimate a quadratic production function with a factor as shifter
germanFarms$decade <- as.factor( c( rep( "70s", 5 ), rep( "80s", 10 ), 
   rep( "90s", 5 ) ) )
estResultShifterFac <- quadFuncEst( yName = "qOutput",
   xNames = c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "decade" ), data = germanFarms )
coef( estResultShifterFac )
print( estResultShifterFac )
# compute fitted values
fitted <- quadFuncCalc( c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "decade" ), data = germanFarms, 
   coef( estResultShifterFac ) )
all.equal( fitted, estResultShifterFac$fitted )

## estimate a quadratic production function with some shifters are logical
estResultShifterLogi2 <- quadFuncEst( yName = "qOutput",
   xNames = c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "time", "reUnif" ), data = germanFarms )
coef( estResultShifterLogi2 )
print( estResultShifterLogi2 )
# compute fitted values
fitted <- quadFuncCalc( c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "time", "reUnif" ), data = germanFarms, 
   coef( estResultShifterLogi2 ) )
all.equal( fitted, estResultShifterLogi2$fitted )

## estimate a quadratic production function with some shifters are factors
estResultShifterFac2 <- quadFuncEst( yName = "qOutput",
   xNames = c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "time", "decade" ), data = germanFarms )
coef( estResultShifterFac2 )
print( estResultShifterFac2 )
# compute fitted values
fitted <- quadFuncCalc( c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "time", "decade" ), data = germanFarms, 
   coef( estResultShifterFac2 ) )
all.equal( fitted, estResultShifterFac2$fitted )

## estimate with further argument passed to lm()
estResult2 <- quadFuncEst( yName = "qOutput",
   xNames = c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, x = TRUE, y = TRUE )
coef( estResult2 )
print( estResult2 )


################ imposing homogeneity #####################
## linear functions
estResultLinHom <- quadFuncEst( yName = "qOutput", 
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   data = germanFarms, linear = TRUE, 
   homWeights = c( qLabor = 0.2, land = 0.5, qVarInput = 0.3 ) )
coef( estResultLinHom )
all.equal( sum( coef( estResultLinHom )[3:4] ), - coef( estResultLinHom )[5],
   check.attributes = FALSE )
vcov( estResultLinHom )
all.equal( rowSums( vcov( estResultLinHom )[ , 3:4 ] ), 
   - vcov( estResultLinHom )[ , 5 ] )
all.equal( colSums( vcov( estResultLinHom )[ 3:4, ] ), 
   - vcov( estResultLinHom )[ 5, ] )
# different order of weights
estResultLinHom2 <- quadFuncEst( yName = "qOutput", 
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   data = germanFarms, linear = TRUE, 
   homWeights = c( qVarInput = 0.3, land = 0.5, qLabor = 0.2 ) )
all.equal( coef( estResultLinHom ), coef( estResultLinHom2 ) )
all.equal( vcov( estResultLinHom ), vcov( estResultLinHom2 ) )
# different order of xNames
estResultLinHom3 <- quadFuncEst( yName = "qOutput", 
   xNames = c( "qLabor", "land", "qVarInput", "time" ),
   data = germanFarms, linear = TRUE, 
   homWeights = c( qLabor = 0.2, land = 0.5, qVarInput = 0.3 ) )
all.equal( coef( estResultLinHom ), 
   coef( estResultLinHom3 )[ c( 1, 5, 2:4, 6:15 ) ],
   check.attributes = FALSE )
all.equal( vcov( estResultLinHom ), 
   vcov( estResultLinHom3 )[ c( 1, 5, 2:4, 6:15 ), c( 1, 5, 2:4, 6:15 ) ],
   check.attributes = FALSE )
# homogenous in all variables
estResultLinHom4 <- quadFuncEst( yName = "qOutput", 
   xNames = c( "time", "qLabor", "land", "qVarInput" ),
   data = germanFarms, linear = TRUE, 
   homWeights = c( qLabor = 0.2, land = 0.5, qVarInput = 0.3, time = 0 ) )
coef( estResultLinHom4 )
all.equal( sum( coef( estResultLinHom4 )[2:4] ), 
 - coef( estResultLinHom4 )[5], check.attributes = FALSE )
vcov( estResultLinHom4 )
all.equal( rowSums( vcov( estResultLinHom4 )[ , 2:4 ] ), 
   - vcov( estResultLinHom4 )[ , 5 ] )
all.equal( colSums( vcov( estResultLinHom4 )[ 2:4, ] ), 
   - vcov( estResultLinHom4 )[ 5, ] )


################ panel data #####################
data( "GrunfeldGreene", package = "systemfit" )
ggData <- plm.data( GrunfeldGreene, c( "firm", "year" ) )
# fixed effects
ggResult <- quadFuncEst( "invest", c( "value", "capital" ), ggData )
coef( ggResult )
print( ggResult )
# random effects
ggResultRan <- quadFuncEst( "invest", c( "value", "capital" ), ggData,
   model = "random", random.method = "amemiya" )
coef( ggResultRan )
print( ggResultRan )

## panel data with a shifter
ggData$yearInt <- as.integer( as.character( ggData$year ) )
ggData$tech <- exp( ggData$yearInt - min( ggData$yearInt ) )
# fixed effects
ggResShifter <- quadFuncEst( "invest", c( "value", "capital" ), ggData,
   shifterNames = "tech" )
coef( ggResShifter )
print.default( ggResShifter )
# random effects
ggResShifterRan <- quadFuncEst( "invest", c( "value", "capital" ), ggData,
   shifterNames = "tech", model = "random", random.method = "amemiya" )
coef( ggResShifterRan )
print.default( ggResShifterRan )

## panel data with a logical variable as shifter
ggData$war <- ggData$yearInt >= 1939 & ggData$yearInt <= 1945
# fixed effects
ggResShifterLogi <- quadFuncEst( "invest", c( "value", "capital" ), ggData,
   shifterNames = "war" )
coef( ggResShifterLogi )
print.default( ggResShifterLogi )
# random effects
ggResShifterLogiRan <- quadFuncEst( "invest", c( "value", "capital" ), ggData,
   shifterNames = "war", model = "random", random.method = "amemiya" )
coef( ggResShifterLogiRan )
print.default( ggResShifterLogiRan )

## panel data with a factor as shifter
ggData$decade <- as.factor( ifelse( ggData$yearInt <= 1939, "30s",
   ifelse( ggData$yearInt <= 1949, "40s", "50s" ) ) )
# fixed effects
ggResShifterFac <- quadFuncEst( "invest", c( "value", "capital" ), ggData,
   shifterNames = "decade" )
coef( ggResShifterFac )
print.default( ggResShifterFac )
# random effects
ggResShifterFacRan <- quadFuncEst( "invest", c( "value", "capital" ), ggData,
   shifterNames = "decade", model = "random", random.method = "amemiya" )
coef( ggResShifterFacRan )
print.default( ggResShifterFacRan )

## linear estimations with panel data
# fixed effects
ggResultLin <- quadFuncEst( "invest", c( "value", "capital" ), ggData,
   linear = TRUE )
coef( ggResultLin )
vcov( ggResultLin )
print( ggResultLin )
# random effects
ggResultLinRan <- quadFuncEst( "invest", c( "value", "capital" ), ggData,
   linear = TRUE, model = "random", random.method = "amemiya" )
coef( ggResultLinRan )
vcov( ggResultLinRan )
print( ggResultLinRan )

## compute partial derivatives of linear estimation results with panel data
# fixed effects
margProducts <- quadFuncDeriv( c( "value", "capital" ),
   data = ggData, coef = coef( ggResultLin ), coefCov = vcov( ggResultLin ) )
sd( margProducts$deriv )
all.equal( margProducts$deriv[1,], coef( ggResultLin )[2:3], 
   check.attributes = FALSE )
all.equal( margProducts$variance[1,], diag( vcov( ggResultLin ) )[2:3], 
   check.attributes = FALSE )
# random effects
margProducts <- quadFuncDeriv( c( "value", "capital" ),
   data = ggData, coef = coef( ggResultLinRan ), coefCov = vcov( ggResultLinRan ) )
sd( margProducts$deriv )
all.equal( margProducts$deriv[1,], coef( ggResultLinRan )[2:3], 
   check.attributes = FALSE )
all.equal( margProducts$variance[1,], diag( vcov( ggResultLinRan ) )[2:3], 
   check.attributes = FALSE )

## imposing homogeneity on linear functions with panel data
ggResultLinHom <- quadFuncEst( "invest", 
   xNames = c( "value", "capital" ), data = ggData,
   linear = TRUE, homWeights = c( value = 0.3, capital = 0.7 ) )
coef( ggResultLinHom )
all.equal( coef( ggResultLinHom )[2], - coef( ggResultLinHom )[3],
   check.attributes = FALSE )
vcov( ggResultLinHom )
all.equal( vcov( ggResultLinHom )[ -1, 2 ], 
   - vcov( ggResultLinHom )[ -1, 3 ] )
all.equal( vcov( ggResultLinHom )[ 2, -1 ], 
   - vcov( ggResultLinHom )[ 3, -1 ] )
