# load the micEconCES package
library( micEconCES )

# seed for random number generation
set.seed( 123 )

# number of observations
nObs <- 200

# create data set with explanatory variables
cesData <- data.frame( xx1 = rchisq( nObs, 10 ), xx2 = rchisq( nObs, 10 ) )

# names of explanatory variables
xxNames <- c( "xx1", "xx2" )

# coefficients
cesCoef <- c( gamma = 1, delta = 0.6, rho = 0.5, phi = 1.1 )

# calculate deterministic endogenous variable
cesData$yd <- cesCalc( xNames = xxNames, data = cesData, coef = cesCoef )
print( cesData$yd )
# check if removing the names of the coefficients makes a difference
all.equal( cesData$yd,
   cesCalc( xNames = xxNames, data = cesData, coef = unname( cesCoef ) ) )
# check if permuting the coefficients makes a difference
all.equal( cesData$yd,
   cesCalc( xNames = xxNames, data = cesData, coef = sample( cesCoef, 4 ) ) )

# adding noise to the endogenous variable
cesData$ys <- cesData$yd + rnorm( nObs )


## Nelder-Mead, CRS
cesNm <- cesEst( "ys", xxNames, cesData )
print.default( cesNm ) 
print( cesNm )
summary( cesNm )
coef( cesNm ) 
vcov( cesNm ) 
coef( summary( cesNm ) )
fitted( cesNm )
residuals( cesNm )

## Nelder-Mead, VRS
cesNmVrs <- cesEst( "ys", xxNames, cesData, vrs = TRUE )
print.default( cesNmVrs )
print( cesNmVrs )
summary( cesNmVrs )
coef( cesNmVrs )
vcov( cesNmVrs )
coef( summary( cesNmVrs ) )
fitted( cesNmVrs )
residuals( cesNmVrs )

## Conjugate Gradients, CRS
cesCg <- cesEst( "ys", xxNames, cesData, method = "CG" )
print.default( cesCg )
print( cesCg )
summary( cesCg )
coef( cesCg )
vcov( cesCg )
coef( summary( cesCg ) )
fitted( cesCg )
residuals( cesCg )

## Conjugate Gradients, VRS
cesCgVrs <- cesEst( "ys", xxNames, cesData, method = "CG", vrs = TRUE )
print.default( cesCgVrs )
print( cesCgVrs )
summary( cesCgVrs )
coef( cesCgVrs )
vcov( cesCgVrs )
coef( summary( cesCgVrs ) )
fitted( cesCgVrs )
residuals( cesCgVrs )

## Simulated Annealing, CRS
cesSann <- cesEst( "ys", xxNames, cesData, method = "SANN" )
print.default( cesSann )
print( cesSann )
summary( cesSann )
coef( cesSann )
vcov( cesSann )
coef( summary( cesSann ) )
fitted( cesSann )
residuals( cesSann )

## Simulated Annealing, VRS
cesSannVrs <- cesEst( "ys", xxNames, cesData, method = "SANN", vrs = TRUE )
print.default( cesSannVrs )
print( cesSannVrs )
summary( cesSannVrs )
coef( cesSannVrs )
vcov( cesSannVrs )
coef( summary( cesSannVrs ) )
fitted( cesSannVrs )
residuals( cesSannVrs )

## BFGS, CRS
cesBfgs <- cesEst( "ys", xxNames, cesData, method = "BFGS" )
print.default( cesBfgs )
print( cesBfgs )
summary( cesBfgs )
coef( cesBfgs )
vcov( cesBfgs )
coef( summary( cesBfgs ) )
fitted( cesBfgs )
residuals( cesBfgs )

## BFGS, VRS
cesBfgsVrs <- cesEst( "ys", xxNames, cesData, method = "BFGS", vrs = TRUE )
print.default( cesBfgsVrs )
print( cesBfgsVrs )
summary( cesBfgsVrs )
coef( cesBfgsVrs )
vcov( cesBfgsVrs )
coef( summary( cesBfgsVrs ) )
fitted( cesBfgsVrs )
residuals( cesBfgsVrs )

## L-BFGS-B with constrained delta, CRS
cesBfgsCon <- cesEst( "ys", xxNames, cesData, method = "L-BFGS-B",
   lower = c( -Inf, 0, -Inf ), upper = c( Inf, 1, Inf ) )
print.default( cesBfgsCon )
print( cesBfgsCon )
summary( cesBfgsCon )
coef( cesBfgsCon )
vcov( cesBfgsCon )
coef( summary( cesBfgsCon ) )
fitted( cesBfgsCon )
residuals( cesBfgsCon )

## L-BFGS-B with constrained delta + phi, VRS
cesBfgsConVrs <- cesEst( "ys", xxNames, cesData, method = "L-BFGS-B",
   vrs = TRUE, lower = c( -Inf, 0, -Inf, 0 ), upper = c( Inf, 1, Inf, Inf ) )
print.default( cesBfgsConVrs )
print( cesBfgsConVrs )
summary( cesBfgsConVrs )
coef( cesBfgsConVrs )
vcov( cesBfgsConVrs )
coef( summary( cesBfgsConVrs ) )
fitted( cesBfgsConVrs )
residuals( cesBfgsConVrs )

## Levenberg-Marquardt, CRS
cesLm <- cesEst( "ys", xxNames, cesData, method = "LM",
   control = nls.lm.control( maxiter = 200 ) )
print.default( cesLm )
print( cesLm )
summary( cesLm )
coef( cesLm )
vcov( cesLm )
coef( summary( cesLm ) )
fitted( cesLm )
residuals( cesLm )

## Levenberg-Marquardt, VRS
cesLmVrs <- cesEst( "ys", xxNames, cesData, method = "LM", vrs = TRUE,
   control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmVrs )
print( cesLmVrs )
summary( cesLmVrs )
coef( cesLmVrs )
vcov( cesLmVrs )
coef( summary( cesLmVrs ) )
fitted( cesLmVrs )
residuals( cesLmVrs )

## Kmenta approximation, CRS
cesKmenta <- cesEst( "ys", xxNames, cesData, method = "Kmenta" )
print.default( cesKmenta )
print( cesKmenta )
summary( cesKmenta )
coef( cesKmenta )
vcov( cesKmenta )
coef( summary( cesKmenta ) )
fitted( cesKmenta )
residuals( cesKmenta )

## Kmenta approximation, VRS
cesKmentaVrs <- cesEst( "ys", xxNames, cesData, method = "Kmenta", vrs = TRUE )
print.default( cesKmentaVrs )
print( cesKmentaVrs )
summary( cesKmentaVrs )
coef( cesKmentaVrs )
vcov( cesKmentaVrs )
coef( summary( cesKmentaVrs ) )
fitted( cesKmentaVrs )
residuals( cesKmentaVrs )
