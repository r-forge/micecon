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

## testing translogEst
# estimate a translog production function
estResult <- translogEst( "qOutput", c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms )

print( estResult )

summary( estResult )

residuals( estResult )

print.default( estResult )

## testing translogCalc
fitted <- translogCalc( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$coef )

all.equal( fitted, estResult$fitted )

## testing translogDeriv
margProducts <- translogDeriv( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$coef, estResult$coefCov )
print( margProducts )

margProductsObs <- translogDeriv( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$coef, estResult$coefCov, yName = "qOutput" )
print( margProductsObs )

germanFarms$fitted <- fitted
margProductsFitted <- translogDeriv( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$coef, estResult$coefCov, yName = "fitted" )
all.equal( margProducts, margProductsFitted )


## testing calculation of elasticities with translogEla
estEla <- translogEla( c( "qLabor", "land", "qVarInput", "time" ), 
   data = germanFarms, coef = coef( estResult ), coefCov = vcov( estResult ) )
print( estEla )


## testing translogHessian
# compute the Hessian matrices
hessians <- translogHessian( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$coef )
print( hessians )

hessiansObs <- translogHessian( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$coef, yName = "qOutput" )
print( hessiansObs )

germanFarms$fitted <- fitted
hessiansFitted <- translogHessian( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$coef, yName = "fitted" )
all.equal( hessians, hessiansFitted )

# compute the bordered Hessian matrices
borderedHessians <- translogHessian( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$coef, bordered = TRUE )
print( borderedHessians )

borderedHessiansObs <- translogHessian( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$coef, yName = "qOutput", bordered = TRUE )
print( borderedHessiansObs )

germanFarms$fitted <- fitted
borderedHessiansFitted <- translogHessian( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$coef, yName = "fitted", bordered = TRUE )
all.equal( borderedHessians, borderedHessiansFitted )

## testing translogCheckMono
test <- translogCheckMono( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ) )
summary( test )
class( test ) <- NULL
print( test )

test <- translogCheckMono( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ), strict = TRUE )
summary( test )
class( test ) <- NULL
print( test )

test <- translogCheckMono( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ), increasing = FALSE )
summary( test )
class( test ) <- NULL
print( test )

test <- translogCheckMono( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ), increasing = FALSE, strict = TRUE )
summary( test )
class( test ) <- NULL
print( test )


## testing translogCheckCurvature
test <- translogCheckCurvature( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ) )
summary( test )
class( test ) <- NULL
print( test )

test <- translogCheckCurvature( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ), quasi = TRUE )
summary( test )
class( test ) <- NULL
print( test )

test <- translogCheckCurvature( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ), convexity = FALSE )
summary( test )
class( test ) <- NULL
print( test )

test <- translogCheckCurvature( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ), convexity = FALSE, quasi = TRUE )
summary( test )
class( test ) <- NULL
print( test )

## testing translogEst with one shifter
germanFarms$tech <- exp( germanFarms$time )
estResultShifter <- translogEst( "qOutput",
   c( "qLabor", "land", "qVarInput" ),
   shifterNames = "tech", data = germanFarms )
print( estResultShifter )
summary( estResultShifter )
residuals( estResultShifter )
print.default( estResultShifter )
# testing translogCalc
fitted <- translogCalc( c( "qLabor", "land", "qVarInput" ),
   shifterNames = "tech", data = germanFarms, coef = estResultShifter$coef )
all.equal( fitted, estResultShifter$fitted )
# testing calculation of elasticities with translogEla
estElaShifter <- translogEla( c( "qLabor", "land", "qVarInput" ), 
   data = germanFarms, coef = coef( estResultShifter ), 
   coefCov = vcov( estResultShifter ) )
print( estElaShifter )

## testing translogEst with two shifters
germanFarms$techSq <- exp( germanFarms$time^2 )
estResultShifter2 <- translogEst( "qOutput",
   c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "tech", "techSq" ), data = germanFarms )
print( estResultShifter2 )
summary( estResultShifter2 )
residuals( estResultShifter2 )
print.default( estResultShifter2 )
# testing translogCalc
fitted <- translogCalc( c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "tech", "techSq" ), data = germanFarms, 
   coef = estResultShifter2$coef )
all.equal( fitted, estResultShifter2$fitted )

## testing translogEst with a logical variable as shifter
germanFarms$reUnif <- germanFarms$time >= 16
estResultShifterLogi <- translogEst( "qOutput",
   c( "qLabor", "land", "qVarInput" ),
   shifterNames = "reUnif", data = germanFarms )
print( estResultShifterLogi )
summary( estResultShifterLogi )
residuals( estResultShifterLogi )
print.default( estResultShifterLogi )
# testing translogCalc
fitted <- translogCalc( c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "reUnif" ), data = germanFarms, 
   coef = estResultShifterLogi$coef )
all.equal( fitted, estResultShifterLogi$fitted )

## testing translogEst with a factor as shifter
germanFarms$decade <- as.factor( c( rep( "70s", 5 ), rep( "80s", 10 ), 
   rep( "90s", 5 ) ) )
estResultShifterFac <- translogEst( "qOutput",
   c( "qLabor", "land", "qVarInput" ),
   shifterNames = "decade", data = germanFarms )
print( estResultShifterFac )
summary( estResultShifterFac )
residuals( estResultShifterFac )
print.default( estResultShifterFac )
# testing translogCalc
fitted <- translogCalc( c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "decade" ), data = germanFarms, 
   coef = estResultShifterFac$coef )
all.equal( fitted, estResultShifterFac$fitted )

## testing translogEst with some shifters and one is logical
estResultShifterLogi2 <- translogEst( "qOutput",
   c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "tech", "reUnif" ), data = germanFarms )
print( estResultShifterLogi2 )
summary( estResultShifterLogi2 )
residuals( estResultShifterLogi2 )
print.default( estResultShifterLogi2 )
# testing translogCalc
fitted <- translogCalc( c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "tech", "reUnif" ), data = germanFarms, 
   coef = estResultShifterLogi2$coef )
all.equal( fitted, estResultShifterLogi2$fitted )

## testing translogEst with some shifters and one is a factor
estResultShifterFac2 <- translogEst( "qOutput",
   c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "tech", "decade" ), data = germanFarms )
print( estResultShifterFac2 )
summary( estResultShifterFac2 )
residuals( estResultShifterFac2 )
print.default( estResultShifterFac2 )
# testing translogCalc
fitted <- translogCalc( c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "tech", "decade" ), data = germanFarms, 
   coef = estResultShifterFac2$coef )
all.equal( fitted, estResultShifterFac2$fitted )

## estimate with further argument passed to lm()
estResult2 <- translogEst( yName = "qOutput",
   xNames = c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, x = TRUE, y = TRUE )
print( estResult2 )
summary( estResult2 )
print.default( estResult2 )

## panel data
data( "GrunfeldGreene", package = "systemfit" )
ggData <- plm.data( GrunfeldGreene, c( "firm", "year" ) )
# fixed effects
ggResult <- translogEst( "invest", c( "value", "capital" ), ggData )
print( ggResult )
print.default( ggResult )
# random effects
ggResultRan <- translogEst( "invest", c( "value", "capital" ), ggData,
   model = "random", random.method = "amemiya" )
print( ggResultRan )
print.default( ggResultRan )

## testing translogEla with panel data
# fixed effects
ggEla <- translogEla( c( "value", "capital" ), 
   data = ggData, coef = coef( ggResult ), 
   coefCov = vcov( ggResult ) )
print( ggEla )
# random effects
ggElaRan <- translogEla( c( "value", "capital" ), 
   data = ggData, coef = coef( ggResultRan ), 
   coefCov = vcov( ggResultRan ) )
print( ggElaRan )

## panel data with a shifter
ggData$yearInt <- as.integer( as.character( ggData$year ) )
ggData$tech <- exp( ggData$yearInt - min( ggData$yearInt ) )
# fixed effects
ggResShifter <- translogEst( "invest", c( "value", "capital" ), ggData,
   shifterNames = "tech" )
print( ggResShifter )
print.default( ggResShifter )
# random effects
ggResShifterRan <- translogEst( "invest", c( "value", "capital" ), ggData,
   shifterNames = "tech", model = "random", random.method = "amemiya" )
print( ggResShifterRan )
print.default( ggResShifterRan )

## panel data with a logical variable as shifter
ggData$war <- ggData$yearInt >= 1939 & ggData$yearInt <= 1945
# fixed effects
ggResShifterLogi <- translogEst( "invest", c( "value", "capital" ), ggData,
   shifterNames = "war" )
print( ggResShifterLogi )
print.default( ggResShifterLogi )
# random effects
ggResShifterLogiRan <- translogEst( "invest", c( "value", "capital" ), ggData,
   shifterNames = "war", model = "random", random.method = "amemiya" )
print( ggResShifterLogiRan )
print.default( ggResShifterLogiRan )

## panel data with a factor as shifter
ggData$decade <- as.factor( ifelse( ggData$yearInt <= 1939, "30s",
   ifelse( ggData$yearInt <= 1949, "40s", "50s" ) ) )
# fixed effects
ggResShifterFac <- translogEst( "invest", c( "value", "capital" ), ggData,
   shifterNames = "decade" )
print( ggResShifterFac )
print.default( ggResShifterFac )
# random effects
ggResShifterFacRan <- translogEst( "invest", c( "value", "capital" ), ggData,
   shifterNames = "decade", model = "random", random.method = "amemiya" )
print( ggResShifterFacRan )
print.default( ggResShifterFacRan )
