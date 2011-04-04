# load the micEconCES package
library( micEconCES )

# seed for random number generation
set.seed( 123 )

# number of observations
nObs <- 200

# create data set with explanatory variables
cesData <- data.frame( xx1 = rchisq( nObs, 10 ), xx2 = rchisq( nObs, 10 ),
   time = c( 0:( nObs - 1 ) ) )

# names of explanatory variables
xxNames <- c( "xx1", "xx2" )

# coefficients
cesCoef <- c( gamma = 1, delta = 0.6, rho = 0.5, nu = 1.1 )
cesCoefTc <- c( cesCoef[ 1 ], lambda = 0.03, cesCoef[ -1 ] )

# calculate deterministic endogenous variable
cesData$yd <- cesCalc( xNames = xxNames, data = cesData, coef = cesCoef )
print( cesData$yd )
# check if removing the names of the coefficients makes a difference
all.equal( cesData$yd,
   cesCalc( xNames = xxNames, data = cesData, coef = unname( cesCoef ) ) )
# check if permuting the coefficients makes a difference
all.equal( cesData$yd,
   cesCalc( xNames = xxNames, data = cesData, coef = sample( cesCoef, 4 ) ) )

# deterministic dependent variable with technological change
cesData$ydTc <- cesCalc( xNames = xxNames, tName = "time", data = cesData, 
   coef = cesCoefTc )
print( cesData$ydTc )
all.equal( cesData$ydTc, 
   cesData$yd * exp( cesCoefTc[ "lambda" ] * c( 0:( nObs - 1 ) ) ) )
# check if removing the names of the coefficients makes a difference
all.equal( cesData$ydTc,
   cesCalc( xNames = xxNames, tName = "time", data = cesData, 
      coef = unname( cesCoefTc ) ) )
# check if permuting the coefficients makes a difference
all.equal( cesData$ydTc,
   cesCalc( xNames = xxNames, tName = "time", data = cesData, 
      coef = cesCoefTc[ c( 4, 1, 3, 5, 2 ) ] ) )


# adding noise to the dependent variables
cesData$ys <- cesData$yd + rnorm( nObs )
print( cesData$ys )

cesData$ysTc <- cesData$ydTc + rnorm( nObs )
print( cesData$ysTc )

cesData$ysMe <- cesData$yd * exp( rnorm( nObs ) )
print( cesData$ysMe )

cesData$ysTcMe <- cesData$ydTc * exp( rnorm( nObs ) )
print( cesData$ysTcMe )


## Nelder-Mead, CRS
cesNm <- cesEst( "ys", xxNames, cesData, method = "Nelder-Mead" )
print.default( cesNm ) 
print( cesNm )
summary( cesNm )
coef( cesNm ) 
vcov( cesNm ) 
coef( summary( cesNm ) )
fitted( cesNm )
residuals( cesNm )

## Nelder-Mead, VRS
cesNmVrs <- cesEst( "ys", xxNames, cesData, vrs = TRUE, method = "NM" )
print.default( cesNmVrs )
print( cesNmVrs )
summary( cesNmVrs )
coef( cesNmVrs )
vcov( cesNmVrs )
coef( summary( cesNmVrs ) )
fitted( cesNmVrs )
residuals( cesNmVrs )

## Nelder-Mead, TC, CRS
cesNmTc <- cesEst( "ysTc", xxNames, tName = "time", data = cesData, 
   method = "Nelder-Mead", control = list( maxit = 1000 ) )
print.default( cesNmTc ) 
print( cesNmTc )
summary( cesNmTc )

## Nelder-Mead, TC, CRS, Start
cesNmTc2 <- cesEst( "ysTc", xxNames, tName = "time", data = cesData, 
   method = "Nelder-Mead", start = c( 2, 0.02, 0.5, 0.25 ) )
print.default( cesNmTc2 ) 
print( cesNmTc2 )
summary( cesNmTc2 )

## Nelder-Mead, TC, VRS
cesNmTcVrs <- cesEst( "ysTc", xxNames, tName = "time", data = cesData, 
   vrs = TRUE, method = "NM", control = list( maxit = 1000 ) )
print.default( cesNmTcVrs )
print( cesNmTcVrs )
summary( cesNmTcVrs )

## Nelder-Mead, TC, VRS, start
cesNmTcVrs2 <- cesEst( "ysTc", xxNames, tName = "time", data = cesData, 
   vrs = TRUE, method = "NM", start = c( 2, 0.02, 0.5, 0.25, 1 ) )
print.default( cesNmTcVrs2 )
print( cesNmTcVrs2 )
summary( cesNmTcVrs2 )

## Nelder-Mead, TC, multErr, VRS
cesNmTcMeVrs <- cesEst( "ysTcMe", xxNames, tName = "time", data = cesData, 
   vrs = TRUE, method = "NM", control = list( maxit = 1000 ), multErr = TRUE )
print.default( cesNmTcMeVrs )
print( cesNmTcMeVrs )
summary( cesNmTcMeVrs )

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

## Conjugate Gradients, CRS, multErr
cesCgMe <- cesEst( "ysMe", xxNames, cesData, method = "CG", multErr = TRUE )
print.default( cesCgMe )
print( cesCgMe )
summary( cesCgMe )
vcov( cesCgMe )

## Simulated Annealing, CRS
# check random number generation
set.seed( 321 )
rnorm( 3 )
cesSann <- cesEst( "ys", xxNames, cesData, method = "SANN" )
print.default( cesSann )
print( cesSann )
summary( cesSann )
coef( cesSann )
vcov( cesSann )
coef( summary( cesSann ) )
fitted( cesSann )
residuals( cesSann )
# check random number generation
rnorm( 4 )

## Simulated Annealing, VRS
# check random number generation
set.seed( 4321 )
rnorm( 4 )
cesSannVrs <- cesEst( "ys", xxNames, cesData, method = "SANN", vrs = TRUE )
print.default( cesSannVrs )
print( cesSannVrs )
summary( cesSannVrs )
coef( cesSannVrs )
vcov( cesSannVrs )
coef( summary( cesSannVrs ) )
fitted( cesSannVrs )
residuals( cesSannVrs )
# check random number generation
rnorm( 5 )

## Simulated Annealing, TC, VRS
cesSannTcVrs <- cesEst( "ysTc", xxNames, tName = "time", data = cesData, 
   method = "SANN", vrs = TRUE )
print.default( cesSannTcVrs )
print( cesSannTcVrs )
summary( cesSannTcVrs )

## Simulated Annealing, TC, multErr
cesSannTcMe <- cesEst( "ysTcMe", xxNames, tName = "time", data = cesData, 
   method = "SANN", multErr = TRUE )
print.default( cesSannTcMe )
print( cesSannTcMe )
summary( cesSannTcMe )
vcov( cesSannTcMe )

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

## BFGS, TC, CRS
cesBfgsTc <- cesEst( "ysTc", xxNames, tName = "time", data = cesData, 
   method = "BFGS" )
print.default( cesBfgsTc )
print( cesBfgsTc )
summary( cesBfgsTc )

## BFGS, TC, VRS
cesBfgsTcVrs <- cesEst( "ysTc", xxNames, tName = "time", data = cesData, 
   method = "BFGS", vrs = TRUE )
print.default( cesBfgsTcVrs )
print( cesBfgsTcVrs )
summary( cesBfgsTcVrs )

## BFGS, multErr, VRS
cesBfgsMeVrs <- cesEst( "ysMe", xxNames, cesData, method = "BFGS", vrs = TRUE,
   multErr = TRUE )
print.default( cesBfgsMeVrs )
print( cesBfgsMeVrs )
summary( cesBfgsMeVrs )
vcov( cesBfgsMeVrs )

## L-BFGS-B with constrained parameters, CRS
cesBfgsCon <- cesEst( "ys", xxNames, cesData, method = "L-BFGS-B" )
print.default( cesBfgsCon )
print( cesBfgsCon )
summary( cesBfgsCon )
coef( cesBfgsCon )
vcov( cesBfgsCon )
coef( summary( cesBfgsCon ) )
fitted( cesBfgsCon )
residuals( cesBfgsCon )

## L-BFGS-B with constrained parameters, VRS
cesBfgsConVrs <- cesEst( "ys", xxNames, cesData, method = "L-BFGS-B",
   vrs = TRUE )
print.default( cesBfgsConVrs )
print( cesBfgsConVrs )
summary( cesBfgsConVrs )
coef( cesBfgsConVrs )
vcov( cesBfgsConVrs )
coef( summary( cesBfgsConVrs ) )
fitted( cesBfgsConVrs )
residuals( cesBfgsConVrs )

## L-BFGS-B with constrained parameters, TC, CRS
cesBfgsConTc <- cesEst( "ysTc", xxNames, tName = "time", data = cesData, 
   method = "L-BFGS-B", control = list( maxit = 500 ) )
print.default( cesBfgsConTc )
print( cesBfgsConTc )
summary( cesBfgsConTc )

## L-BFGS-B with constrained parameters, TC, VRS
cesBfgsConTcVrs <- cesEst( "ysTc", xxNames, tName = "time", data = cesData, 
   method = "L-BFGS-B", vrs = TRUE )
print.default( cesBfgsConTcVrs )
print( cesBfgsConTcVrs )
summary( cesBfgsConTcVrs )

## L-BFGS-B with constrained parameters, TC, multErr, VRS
cesBfgsConTcMeVrs <- cesEst( "ysTcMe", xxNames, tName = "time", data = cesData, 
   method = "L-BFGS-B", vrs = TRUE, multErr = TRUE )
print.default( cesBfgsConTcMeVrs )
print( cesBfgsConTcMeVrs )
summary( cesBfgsConTcMeVrs )

## Levenberg-Marquardt, CRS
cesLm <- cesEst( "ys", xxNames, cesData,
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
cesLmVrs <- cesEst( "ys", xxNames, cesData, vrs = TRUE,
   control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmVrs )
print( cesLmVrs )
summary( cesLmVrs )
coef( cesLmVrs )
vcov( cesLmVrs )
coef( summary( cesLmVrs ) )
fitted( cesLmVrs )
residuals( cesLmVrs )

## Levenberg-Marquardt, TC, CRS
cesLmTc <- cesEst( "ysTc", xxNames, tName = "time", data = cesData )
print.default( cesLmTc )
print( cesLmTc )
summary( cesLmTc )

## Levenberg-Marquardt, TC, VRS
cesLmTcVrs <- cesEst( "ysTc", xxNames, tName = "time", data = cesData, 
   vrs = TRUE )
print.default( cesLmTcVrs )
print( cesLmTcVrs )
summary( cesLmTcVrs )

## Levenberg-Marquardt, CRS, multErr
cesLmMe <- cesEst( "ysMe", xxNames, cesData, multErr = TRUE,
   control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmMe )
print( cesLmMe )
summary( cesLmMe )
vcov( cesLmMe )

## Newton-type, CRS
cesNewton <- cesEst( "ys", xxNames, cesData, method = "Newton" )
print.default( cesNewton )
print( cesNewton )
summary( cesNewton )
coef( cesNewton )
vcov( cesNewton )
coef( summary( cesNewton ) )
fitted( cesNewton )
residuals( cesNewton )

## Newton-type, VRS
cesNewtonVrs <- cesEst( "ys", xxNames, cesData, method = "Newton", vrs = TRUE )
print.default( cesNewtonVrs )
print( cesNewtonVrs )
summary( cesNewtonVrs )
coef( cesNewtonVrs )
vcov( cesNewtonVrs )
coef( summary( cesNewtonVrs ) )
fitted( cesNewtonVrs )
residuals( cesNewtonVrs )

## Newton-type, TC, CRS
cesNewtonTc <- cesEst( "ysTc", xxNames, tName = "time", data = cesData, 
   method = "Newton" )
print.default( cesNewtonTc )
print( cesNewtonTc )
summary( cesNewtonTc )

## Newton-type, TC, VRS
cesNewtonTcVrs <- cesEst( "ysTc", xxNames, tName = "time", data = cesData, 
   method = "Newton", vrs = TRUE, iterlim = 500 )
print.default( cesNewtonTcVrs )
print( cesNewtonTcVrs )
summary( cesNewtonTcVrs )

## Newton-type, TC, multErr, VRS
cesNewtonTcMeVrs <- cesEst( "ysTcMe", xxNames, tName = "time", data = cesData, 
   method = "Newton", vrs = TRUE, multErr = TRUE, iterlim = 500 )
print.default( cesNewtonTcMeVrs )
print( cesNewtonTcMeVrs )
summary( cesNewtonTcMeVrs )
vcov( cesNewtonTcMeVrs )

## PORT, CRS
cesPort <- cesEst( "ys", xxNames, cesData, method = "PORT" )
print.default( cesPort )
print( cesPort )
summary( cesPort )
coef( cesPort )
vcov( cesPort )
coef( summary( cesPort ) )
fitted( cesPort )
residuals( cesPort )

## PORT, VRS
cesPortVrs <- cesEst( "ys", xxNames, cesData, method = "PORT", vrs = TRUE )
print.default( cesPortVrs )
print( cesPortVrs )
summary( cesPortVrs )
coef( cesPortVrs )
vcov( cesPortVrs )
coef( summary( cesPortVrs ) )
fitted( cesPortVrs )
residuals( cesPortVrs )

## PORT, TC, CRS
cesPortTc <- cesEst( "ysTc", xxNames, tName = "time", data = cesData, 
   method = "PORT" )
print.default( cesPortTc )
print( cesPortTc )
summary( cesPortTc )

## PORT, TC, VRS
cesPortTcVrs <- cesEst( "ysTc", xxNames, tName = "time", data = cesData, 
   method = "PORT", vrs = TRUE )
print.default( cesPortTcVrs )
print( cesPortTcVrs )
summary( cesPortTcVrs )

## PORT, multErr, VRS
cesPortMeVrs <- cesEst( "ysMe", xxNames, cesData, method = "PORT", vrs = TRUE,
   multErr = TRUE )
print.default( cesPortMeVrs )
print( cesPortMeVrs )
summary( cesPortMeVrs )
vcov( cesPortMeVrs )

## DE, CRS
# check random number generation
set.seed( 54321 )
rnorm( 5 )
cesDe <- cesEst( "ys", xxNames, cesData, method = "DE",
   control = DEoptim.control( trace = FALSE ) )
print.default( cesDe )
print( cesDe )
summary( cesDe )
coef( cesDe )
vcov( cesDe )
coef( summary( cesDe ) )
fitted( cesDe )
residuals( cesDe )
print( fitted( cesDe ) + residuals( cesDe ) )
# check random number generation
rnorm( 4 )

## DE, VRS
# check random number generation
set.seed( 654321 )
rnorm( 4 )
cesDeVrs <- cesEst( "ys", xxNames, cesData, method = "DE", vrs = TRUE,
   control = DEoptim.control( trace = FALSE ) )
print.default( cesDeVrs )
print( cesDeVrs )
summary( cesDeVrs )
coef( cesDeVrs )
vcov( cesDeVrs )
coef( summary( cesDeVrs ) )
fitted( cesDeVrs )
residuals( cesDeVrs )
print( fitted( cesDeVrs ) + residuals( cesDeVrs ) )
# check random number generation
rnorm( 5 )

## DE, TC, VRS
cesDeTcVrs <- cesEst( "ysTc", xxNames, tName = "time", data = cesData, 
   method = "DE", vrs = TRUE, control = DEoptim.control( trace = FALSE ) )
print.default( cesDeTcVrs )
print( cesDeTcVrs )
summary( cesDeTcVrs )

## DE, TC, multErr
cesDeTcMe <- cesEst( "ysMe", xxNames, tName = "time", data = cesData, 
   method = "DE", multErr = TRUE, control = DEoptim.control( trace = FALSE ) )
print.default( cesDeTcMe )
print( cesDeTcMe )
summary( cesDeTcMe )

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

## nls, CRS
cesNls <- cesEst( "ys", xxNames, cesData, method = "nls" )
print.default( cesNls )
print( cesNls )
summary( cesNls )
coef( cesNls )
vcov( cesNls )
coef( summary( cesNls ) )
fitted( cesNls )
residuals( cesNls )

## nls, VRS
cesNlsVrs <- cesEst( "ys", xxNames, cesData, method = "nls", vrs = TRUE )
print.default( cesNlsVrs )
print( cesNlsVrs )
summary( cesNlsVrs )
coef( cesNlsVrs )
vcov( cesNlsVrs )
coef( summary( cesNlsVrs ) )
fitted( cesNlsVrs )
residuals( cesNlsVrs )

## nls, TC, VRS
cesNlsTcVrs <- cesEst( "ysTc", xxNames, tName = "time", data = cesData, 
   method = "nls", vrs = TRUE )
print.default( cesNlsTcVrs )
print( cesNlsTcVrs )
summary( cesNlsTcVrs )

## nls, multErr, VRS
cesNlsMeVrs <- cesEst( "ysMe", xxNames, cesData, method = "nls", vrs = TRUE,
   multErr = TRUE )
print.default( cesNlsMeVrs )
print( cesNlsMeVrs )
summary( cesNlsMeVrs )

## nls, TC, multErr
cesNlsTcMe <- cesEst( "ysTcMe", xxNames, tName = "time", data = cesData, 
   method = "nls", multErr = TRUE )
print.default( cesNlsTcMe )
print( cesNlsTcMe )
summary( cesNlsTcMe )

## nls, TC, multErr, VRS
cesNlsTcMeVrs <- try( cesEst( "ysTcMe", xxNames, tName = "time", data = cesData, 
   method = "nls", vrs = TRUE, multErr = TRUE ) )
print( cesNlsTcMeVrs )

