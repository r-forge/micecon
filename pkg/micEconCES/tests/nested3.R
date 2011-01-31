library( "micEconCES" )
library( "maxLik" )

data( "MishraCES" )

bOld <- c( "gamma_1" = 2, "gamma_2" = 200, "delta_1" = 0.6, "delta_2" = 0.7, 
   "rho_1" = 0.5, "rho" = 0.6 )
bOldVrs <- c( bOld, "nu" = 1.1 )
xNames <- c( "X1", "X2", "X3" )

# normalize gamma_1 to 1
b <- bOld
b[ "gamma_2" ] <- bOld[ "gamma_2" ] * 
   ( bOld[ "delta_2" ] * bOld[ "gamma_1" ]^( - bOld[ "rho" ] ) + ( 1 - bOld[ "delta_2" ] ) )^( 
      -1 / bOld[ "rho" ] )
b[ "delta_2" ] <- bOld[ "delta_2" ] * bOld[ "gamma_1" ]^( - bOld[ "rho" ] ) /
   ( bOld[ "delta_2" ] * bOld[ "gamma_1" ]^( - bOld[ "rho" ] ) + ( 1 - bOld[ "delta_2" ] ) )
b <- b[ names( b ) != "gamma_1" ]
rm( bOld )

bVrs <- bOldVrs
bVrs[ "gamma_2" ] <- bOldVrs[ "gamma_2" ] * 
   ( bOldVrs[ "delta_2" ] * bOldVrs[ "gamma_1" ]^( - bOldVrs[ "rho" ] ) + ( 1 - bOldVrs[ "delta_2" ] ) )^( 
      - bOldVrs[ "nu" ] / bOldVrs[ "rho" ] )
bVrs[ "delta_2" ] <- bOldVrs[ "delta_2" ] * bOldVrs[ "gamma_1" ]^( - bOldVrs[ "rho" ] ) /
   ( bOldVrs[ "delta_2" ] * bOldVrs[ "gamma_1" ]^( - bOldVrs[ "rho" ] ) + ( 1 - bOldVrs[ "delta_2" ] ) )
bVrs <- bVrs[ names( bVrs ) != "gamma_1" ]
rm( bOldVrs )


######################### checking cesCalc() ###############################
MishraCES$Y3 <- cesCalc( xNames = xNames, 
   data = MishraCES, coef = b, nested = TRUE )
MishraCES$Y3

# VRS
MishraCES$Y3Vrs <- cesCalc( xNames = xNames, 
   data = MishraCES, coef = bVrs, nested = TRUE )
MishraCES$Y3Vrs


## check cesCalc() with rho equal to zero and close to zero
# vector with values around 0 for coefficient "rho"
rhos <- c( -exp(-(1:20)),0,exp(-(20:1)) )
# rhos <- c( -2^(-(1:40)),0,2^(-(40:1)) )

# matrix for returned endogenous variables
yRho <- matrix( NA, nrow = length( rhos ), ncol = 7 )
rownames( yRho ) <- c( -(1:20), 0, (20:1) )

# calculate endogenous variables
coefRho <- b
for( i in 1:length( rhos ) ) {
   coefRho[ "rho" ] <- rhos[ i ]
   yRho[ i, ] <- cesCalc( xNames = xNames, data = MishraCES[1:7,], coef = coefRho,
      nested = TRUE )
}

# print "raw" endogenous values
print( yRho )

# print endogenous variables for different rhos (adjusted with the y at rho=0)
for( i in 1:ncol( yRho ) ) {
   print( format( round( yRho[ , i, drop = FALSE ] - yRho[21,i], 11 ),
      scientific = FALSE ) )
}

## check cesCalc() with rho_1 equal to zero and close to zero
# matrix for returned endogenous variables
yRho1 <- matrix( NA, nrow = length( rhos ), ncol = 7 )
rownames( yRho1 ) <- c( -(1:20), 0, (20:1) )

# calculate endogenous variables
coefRho1 <- b
for( i in 1:length( rhos ) ) {
   coefRho1[ "rho_1" ] <- rhos[ i ]
   yRho1[ i, ] <- cesCalc( xNames = xNames, data = MishraCES[1:7,], coef = coefRho1,
      nested = TRUE )
}

# print "raw" endogenous values
print( yRho1 )


## check cesCalc() with rho_1 and rho equal to zero and close to zero
# array for returned endogenous variables
yRho2 <- array( NA, c( length( rhos ), length( rhos ), 2 ) )
dimnames( yRho2 ) <- list( -20:20, -20:20, 1:2 )

# calculate endogenous variables
coefRho2 <- bVrs
for( i in 1:length( rhos ) ) {
   for( j in 1:length( rhos ) ) {
      coefRho2[ "rho" ] <- rhos[ i ]
      coefRho2[ "rho_1" ] <- rhos[ j ]
      yRho2[ i, j, ] <- cesCalc( xNames = xNames, data = MishraCES[ 1:2, ], 
         coef = coefRho2, nested = TRUE )
   }
}

# print "raw" endogenous values
print( yRho2 )


########################## checking cesDerivCoef ##############################
cesDeriv <- micEconCES:::cesDerivCoef( par = b, xNames = xNames, 
   data = MishraCES, vrs = FALSE, nested = TRUE, 
   rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
f <- function( par ) {
   return( cesCalc( xNames = xNames, data = MishraCES, coef = par, 
      nested = TRUE ) )
}
cesDerivNum <- numericGradient( f, t0 = b )
all.equal( cesDeriv, cesDerivNum )
print( cesDeriv )

# VRS
cesDerivVrs <- micEconCES:::cesDerivCoef( par = bVrs, xNames = xNames, 
   data = MishraCES, vrs = TRUE, nested = TRUE,
   rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
cesDerivVrsNum <- numericGradient( f, t0 = bVrs )
all.equal( cesDerivVrs, cesDerivVrsNum )
print( cesDerivVrs )


## check cesDerivCoef() with rho equal to zero and close to zero
# array for returned partial derivatives
deriv <- array( NA, c( length( rhos ), 5, length( bVrs ) ) )
dimnames( deriv ) <- list( rhos, 1:5, names( bVrs ) )

coefRho <- bVrs
# calculate the derivatives
for( i in 1:length( rhos ) ) {
   # coefficients
   coefRho[ "rho" ] <- rhos[ i ]
   deriv[ i, , ] <- micEconCES:::cesDerivCoef( par = coefRho,
      xNames = xNames, data = MishraCES[1:5,], nested = TRUE, vrs = TRUE,
      rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
}

# print array of derivatives
print( deriv )

# print derivatives for different rhos (adjusted with the derivatives at rho=0)
for( k in 1:dim( deriv )[3] ) {
   for( i in 1:5 ) {
      print( format( round( deriv[ , i, k, drop = FALSE ] -
       deriv[ 21, i, k ], 11 ), scientific = FALSE ) )
   }
}

## check cesDerivCoef() with rho_1 equal to zero and close to zero
# array for returned partial derivatives
deriv1 <- array( NA, c( length( rhos ), 5, length( bVrs ) ) )
dimnames( deriv1 ) <- list( rhos, 1:5, names( bVrs ) )

coefRho1 <- bVrs
# calculate the derivatives
for( i in 1:length( rhos ) ) {
   # coefficients
   coefRho1[ "rho_1" ] <- rhos[ i ]
   deriv1[ i, , ] <- micEconCES:::cesDerivCoef( par = coefRho1,
      xNames = xNames, data = MishraCES[1:5,], nested = TRUE, vrs = TRUE,
      rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
}

# print array of derivatives
print( deriv1 )


## check cesDerivCoef() with rho and rho_1 equal to zero and close to zero
# array for returned partial derivatives
deriv2 <- array( NA, c( length( rhos ), length( rhos ), 2, length( bVrs ) ) )
dimnames( deriv2 ) <- list( c(-20:20), c(-20:20), 1:2, names( bVrs ) )

coefRho2 <- bVrs
# calculate the derivatives
for( i in 1:length( rhos ) ) {
   for( j in 1:length( rhos ) ) {
      # coefficients
      coefRho2[ "rho" ] <- rhos[ i ]
      coefRho2[ "rho_1" ] <- rhos[ j ]
      deriv2[ i, j, , ] <- micEconCES:::cesDerivCoef( par = coefRho2,
         xNames = xNames, data = MishraCES[1:2,], nested = TRUE, vrs = TRUE,
         rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
   }
}

# print array of derivatives
print( deriv2 )



########################## checking cesEst ###################################
set.seed( 345 )
MishraCES$yObs <- MishraCES$Y3 + 400 * rnorm( nrow( MishraCES ) )

## Nelder-Mead, CRS
cesNm <- cesEst( "yObs", xNames, data = MishraCES, method = "Nelder-Mead" )
print.default( cesNm ) 
print( cesNm )
summary( cesNm )
coef( cesNm ) 
vcov( cesNm ) 
coef( summary( cesNm ) )
fitted( cesNm )
residuals( cesNm )

## Nelder-Mead, VRS
cesNmVrs <- cesEst( "yObs", xNames, data = MishraCES, vrs = TRUE, method = "NM",
   control = list( maxit = 1000 ) )
print.default( cesNmVrs )
print( cesNmVrs )
summary( cesNmVrs )
coef( cesNmVrs )
vcov( cesNmVrs )
coef( summary( cesNmVrs ) )
fitted( cesNmVrs )
residuals( cesNmVrs )

## Conjugate Gradients, CRS
cesCg <- cesEst( "yObs", xNames, data = MishraCES, method = "CG" )
print.default( cesCg )
print( cesCg )
summary( cesCg )
coef( cesCg )
vcov( cesCg )
coef( summary( cesCg ) )
fitted( cesCg )
residuals( cesCg )

## Conjugate Gradients, VRS
cesCgVrs <- cesEst( "yObs", xNames, data = MishraCES, method = "CG", 
   vrs = TRUE )
print.default( cesCgVrs )
print( cesCgVrs )
summary( cesCgVrs )
coef( cesCgVrs )
vcov( cesCgVrs )
coef( summary( cesCgVrs ) )
fitted( cesCgVrs )
residuals( cesCgVrs )

## Simulated Annealing, CRS
cesSann <- cesEst( "yObs", xNames, data = MishraCES, method = "SANN" )
print.default( cesSann )
print( cesSann )
summary( cesSann )
coef( cesSann )
vcov( cesSann )
coef( summary( cesSann ) )
fitted( cesSann )
residuals( cesSann )

## Simulated Annealing, VRS
cesSannVrs <- cesEst( "yObs", xNames, data = MishraCES, method = "SANN", 
   vrs = TRUE )
print.default( cesSannVrs )
print( cesSannVrs )
summary( cesSannVrs )
coef( cesSannVrs )
vcov( cesSannVrs )
coef( summary( cesSannVrs ) )
fitted( cesSannVrs )
residuals( cesSannVrs )

## BFGS, CRS
cesBfgs <- cesEst( "yObs", xNames, data = MishraCES, method = "BFGS",
   control = list( maxit = 500 ) )
print.default( cesBfgs )
print( cesBfgs )
summary( cesBfgs )
coef( cesBfgs )
vcov( cesBfgs )
coef( summary( cesBfgs ) )
fitted( cesBfgs )
residuals( cesBfgs )

## BFGS, VRS
cesBfgsVrs <- cesEst( "yObs", xNames, data = MishraCES, method = "BFGS", 
   vrs = TRUE, control = list( maxit = 500 ) )
print.default( cesBfgsVrs )
print( cesBfgsVrs )
summary( cesBfgsVrs )
coef( cesBfgsVrs )
vcov( cesBfgsVrs )
coef( summary( cesBfgsVrs ) )
fitted( cesBfgsVrs )
residuals( cesBfgsVrs )

## L-BFGS-B with constrained parameters, CRS
cesBfgsCon <- cesEst( "yObs", xNames, data = MishraCES, method = "L-BFGS-B" )
print.default( cesBfgsCon )
print( cesBfgsCon )
summary( cesBfgsCon )
coef( cesBfgsCon )
vcov( cesBfgsCon )
coef( summary( cesBfgsCon ) )
fitted( cesBfgsCon )
residuals( cesBfgsCon )

## L-BFGS-B with constrained parameters, VRS
cesBfgsConVrs <- cesEst( "yObs", xNames, data = MishraCES, method = "L-BFGS-B",
   vrs = TRUE, control = list( maxit = 500 ) )
print.default( cesBfgsConVrs )
print( cesBfgsConVrs )
summary( cesBfgsConVrs )
coef( cesBfgsConVrs )
vcov( cesBfgsConVrs )
coef( summary( cesBfgsConVrs ) )
fitted( cesBfgsConVrs )
residuals( cesBfgsConVrs )

## Levenberg-Marquardt, CRS
cesLm <- cesEst( "yObs", xNames, data = MishraCES,
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
cesLmVrs <- cesEst( "yObs", xNames, data = MishraCES, vrs = TRUE,
   control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmVrs )
print( cesLmVrs )
summary( cesLmVrs )
coef( cesLmVrs )
vcov( cesLmVrs )
coef( summary( cesLmVrs ) )
fitted( cesLmVrs )
residuals( cesLmVrs )

## Newton-type, CRS
cesNewton <- cesEst( "yObs", xNames, data = MishraCES, method = "Newton" )
print.default( cesNewton )
print( cesNewton )
summary( cesNewton )
coef( cesNewton )
vcov( cesNewton )
coef( summary( cesNewton ) )
fitted( cesNewton )
residuals( cesNewton )

## Newton-type, VRS
cesNewtonVrs <- cesEst( "yObs", xNames, data = MishraCES, method = "Newton", 
   vrs = TRUE )
print.default( cesNewtonVrs )
print( cesNewtonVrs )
summary( cesNewtonVrs )
coef( cesNewtonVrs )
vcov( cesNewtonVrs )
coef( summary( cesNewtonVrs ) )
fitted( cesNewtonVrs )
residuals( cesNewtonVrs )

## PORT, CRS
cesPort <- cesEst( "yObs", xNames, data = MishraCES, method = "PORT",
   control = list( eval.max = 500 ) )
print.default( cesPort )
print( cesPort )
summary( cesPort )
coef( cesPort )
vcov( cesPort )
coef( summary( cesPort ) )
fitted( cesPort )
residuals( cesPort )

## PORT, VRS
cesPortVrs <- cesEst( "yObs", xNames, data = MishraCES, method = "PORT", 
   vrs = TRUE, control = list( eval.max = 500, iter.max = 500 ) )
print.default( cesPortVrs )
print( cesPortVrs )
summary( cesPortVrs )
coef( cesPortVrs )
vcov( cesPortVrs )
coef( summary( cesPortVrs ) )
fitted( cesPortVrs )
residuals( cesPortVrs )

## DE, CRS
cesDe <- cesEst( "yObs", xNames, data = MishraCES, method = "DE",
   control = DEoptim.control( trace = FALSE, NP = 60 ) )
print.default( cesDe )
print( cesDe )
summary( cesDe )
coef( cesDe )
vcov( cesDe )
coef( summary( cesDe ) )
fitted( cesDe )
residuals( cesDe )
print( fitted( cesDe ) + residuals( cesDe ) )

## DE, VRS
cesDeVrs <- cesEst( "yObs", xNames, data = MishraCES, method = "DE", vrs = TRUE,
   control = DEoptim.control( trace = FALSE, NP = 70 ) )
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


########## Fixing Rho_1 ##############
## Levenberg-Marquardt, rho_1 fixed, CRS
cesLmRho1 <- cesEst( "yObs", xNames, data = MishraCES, rho1 = 1,
   control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmRho1 )
print( cesLmRho1 )
summary( cesLmRho1 )
coef( cesLmRho1 )
vcov( cesLmRho1 )
coef( summary( cesLmRho1 ) )
fitted( cesLmRho1 )
residuals( cesLmRho1 )

## Levenberg-Marquardt, rho_1 fixed, VRS
cesLmVrsRho1 <- cesEst( "yObs", xNames, data = MishraCES, vrs = TRUE,
   rho1 = 0, control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmVrsRho1 )
print( cesLmVrsRho1 )
summary( cesLmVrsRho1 )
coef( cesLmVrsRho1 )
vcov( cesLmVrsRho1 )
coef( summary( cesLmVrsRho1 ) )
fitted( cesLmVrsRho1 )
residuals( cesLmVrsRho1 )

## PORT, rho_1 fixed, CRS
cesPortRho1 <- cesEst( "yObs", xNames, data = MishraCES, method = "PORT",
   rho1 = 0, control = list( eval.max = 500 ) )
print.default( cesPortRho1 )
print( cesPortRho1 )
summary( cesPortRho1 )
coef( cesPortRho1 )
vcov( cesPortRho1 )
coef( summary( cesPortRho1 ) )
fitted( cesPortRho1 )
residuals( cesPortRho1 )

## PORT, rho_1 fixed, VRS
cesPortVrsRho1 <- cesEst( "yObs", xNames, data = MishraCES, method = "PORT", 
   vrs = TRUE, rho1 = 1, 
   control = list( eval.max = 500, iter.max = 500 ) )
print.default( cesPortVrsRho1 )
print( cesPortVrsRho1 )
summary( cesPortVrsRho1 )
coef( cesPortVrsRho1 )
vcov( cesPortVrsRho1 )
coef( summary( cesPortVrsRho1 ) )
fitted( cesPortVrsRho1 )
residuals( cesPortVrsRho1 )


########## Fixing Rho ##############
## Levenberg-Marquardt, rho fixed, CRS
cesLmRho <- cesEst( "yObs", xNames, data = MishraCES, rho = 1,
   control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmRho )
print( cesLmRho )
summary( cesLmRho )
coef( cesLmRho )
vcov( cesLmRho )
coef( summary( cesLmRho ) )
fitted( cesLmRho )
residuals( cesLmRho )

## Levenberg-Marquardt, rho fixed, VRS
cesLmVrsRho <- cesEst( "yObs", xNames, data = MishraCES, vrs = TRUE,
   rho = 0, control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmVrsRho )
print( cesLmVrsRho )
summary( cesLmVrsRho )
coef( cesLmVrsRho )
vcov( cesLmVrsRho )
coef( summary( cesLmVrsRho ) )
fitted( cesLmVrsRho )
residuals( cesLmVrsRho )

## PORT, rho fixed, CRS
cesPortRho <- cesEst( "yObs", xNames, data = MishraCES, method = "PORT",
   rho = 0, control = list( eval.max = 500 ) )
print.default( cesPortRho )
print( cesPortRho )
summary( cesPortRho )
coef( cesPortRho )
vcov( cesPortRho )
coef( summary( cesPortRho ) )
fitted( cesPortRho )
residuals( cesPortRho )

## PORT, rho fixed, VRS
cesPortVrsRho <- cesEst( "yObs", xNames, data = MishraCES, method = "PORT", 
   vrs = TRUE, rho = 1, 
   control = list( eval.max = 500, iter.max = 500 ) )
print.default( cesPortVrsRho )
print( cesPortVrsRho )
summary( cesPortVrsRho )
coef( cesPortVrsRho )
vcov( cesPortVrsRho )
coef( summary( cesPortVrsRho ) )
fitted( cesPortVrsRho )
residuals( cesPortVrsRho )


########## Fixing Rho_1 and Rho ##############
## Levenberg-Marquardt, rho_1 and rho fixed, CRS
cesLmRho2 <- cesEst( "yObs", xNames, data = MishraCES, rho1 = 1, rho = 1,
   control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmRho2 )
print( cesLmRho2 )
summary( cesLmRho2 )
coef( cesLmRho2 )
vcov( cesLmRho2 )
coef( summary( cesLmRho2 ) )
fitted( cesLmRho2 )
residuals( cesLmRho2 )

## Levenberg-Marquardt, rho_1 and rho fixed, VRS
cesLmVrsRho2 <- cesEst( "yObs", xNames, data = MishraCES, vrs = TRUE,
   rho1 = 0, rho = 0, control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmVrsRho2 )
print( cesLmVrsRho2 )
summary( cesLmVrsRho2 )
coef( cesLmVrsRho2 )
vcov( cesLmVrsRho2 )
coef( summary( cesLmVrsRho2 ) )
fitted( cesLmVrsRho2 )
residuals( cesLmVrsRho2 )

## PORT, rho_1 and rho fixed, CRS
cesPortRho2 <- cesEst( "yObs", xNames, data = MishraCES, method = "PORT",
   rho1 = 0, rho = 1, control = list( eval.max = 500 ) )
print.default( cesPortRho2 )
print( cesPortRho2 )
summary( cesPortRho2 )
coef( cesPortRho2 )
vcov( cesPortRho2 )
coef( summary( cesPortRho2 ) )
fitted( cesPortRho2 )
residuals( cesPortRho2 )

## PORT, rho_1 and rho fixed, VRS
cesPortVrsRho2 <- cesEst( "yObs", xNames, data = MishraCES, method = "PORT", 
   vrs = TRUE, rho1 = 1, rho = 0,
   control = list( eval.max = 500, iter.max = 500 ) )
print.default( cesPortVrsRho2 )
print( cesPortVrsRho2 )
summary( cesPortVrsRho2 )
coef( cesPortVrsRho2 )
vcov( cesPortVrsRho2 )
coef( summary( cesPortVrsRho2 ) )
fitted( cesPortVrsRho2 )
residuals( cesPortVrsRho2 )


########## Grid Search for Rho_1 ##############
## Levenberg-Marquardt, line search for rho_1, CRS
gridRhos <- (1:12)/6-0.4
gridRhos <- c( gridRhos[ 1:2 ], 0, gridRhos[ -c(1:2) ] )
cesLmGrid1 <- cesEst( "yObs", xNames, data = MishraCES, rho1 = gridRhos,
   control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmGrid1 )
print( cesLmGrid1 )
summary( cesLmGrid1 )
coef( cesLmGrid1 )
vcov( cesLmGrid1 )
coef( summary( cesLmGrid1 ) )
fitted( cesLmGrid1 )
residuals( cesLmGrid1 )

## Levenberg-Marquardt, line search for rho_1, VRS
cesLmVrsGrid1 <- cesEst( "yObs", xNames, data = MishraCES, vrs = TRUE,
   rho1 = gridRhos, control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmVrsGrid1 )
print( cesLmVrsGrid1 )
summary( cesLmVrsGrid1 )
coef( cesLmVrsGrid1 )
vcov( cesLmVrsGrid1 )
coef( summary( cesLmVrsGrid1 ) )
fitted( cesLmVrsGrid1 )
residuals( cesLmVrsGrid1 )

## PORT, line search for rho_1, CRS
cesPortGrid1 <- cesEst( "yObs", xNames, data = MishraCES, method = "PORT",
   rho1 = gridRhos, control = list( eval.max = 500 ) )
print.default( cesPortGrid1 )
print( cesPortGrid1 )
summary( cesPortGrid1 )
coef( cesPortGrid1 )
vcov( cesPortGrid1 )
coef( summary( cesPortGrid1 ) )
fitted( cesPortGrid1 )
residuals( cesPortGrid1 )

## PORT, line search for rho_1, VRS
cesPortVrsGrid1 <- cesEst( "yObs", xNames, data = MishraCES, method = "PORT", 
   vrs = TRUE, rho1 = gridRhos, 
   control = list( eval.max = 500, iter.max = 500 ) )
print.default( cesPortVrsGrid1 )
print( cesPortVrsGrid1 )
summary( cesPortVrsGrid1 )
coef( cesPortVrsGrid1 )
vcov( cesPortVrsGrid1 )
coef( summary( cesPortVrsGrid1 ) )
fitted( cesPortVrsGrid1 )
residuals( cesPortVrsGrid1 )


########## Grid Search for Rho ##############
## Levenberg-Marquardt, Grid Search, CRS
cesLmGrid <- cesEst( "yObs", xNames, data = MishraCES, rho = gridRhos,
   control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmGrid )
print( cesLmGrid )
summary( cesLmGrid )
coef( cesLmGrid )
vcov( cesLmGrid )
coef( summary( cesLmGrid ) )
fitted( cesLmGrid )
residuals( cesLmGrid )

## Levenberg-Marquardt, Grid Search, VRS
cesLmVrsGrid <- cesEst( "yObs", xNames, data = MishraCES, vrs = TRUE,
   rho = gridRhos, control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmVrsGrid )
print( cesLmVrsGrid )
summary( cesLmVrsGrid )
coef( cesLmVrsGrid )
vcov( cesLmVrsGrid )
coef( summary( cesLmVrsGrid ) )
fitted( cesLmVrsGrid )
residuals( cesLmVrsGrid )

## PORT, Grid Search, CRS
cesPortGrid <- cesEst( "yObs", xNames, data = MishraCES, method = "PORT",
   rho = gridRhos, control = list( eval.max = 500 ) )
print.default( cesPortGrid )
print( cesPortGrid )
summary( cesPortGrid )
coef( cesPortGrid )
vcov( cesPortGrid )
coef( summary( cesPortGrid ) )
fitted( cesPortGrid )
residuals( cesPortGrid )

## PORT, Grid Search, VRS
cesPortVrsGrid <- cesEst( "yObs", xNames, data = MishraCES, method = "PORT", 
   vrs = TRUE, rho = gridRhos, 
   control = list( eval.max = 500, iter.max = 500 ) )
print.default( cesPortVrsGrid )
print( cesPortVrsGrid )
summary( cesPortVrsGrid )
coef( cesPortVrsGrid )
vcov( cesPortVrsGrid )
coef( summary( cesPortVrsGrid ) )
fitted( cesPortVrsGrid )
residuals( cesPortVrsGrid )


########## Grid Search for rho_1 and rho ##############
## Levenberg-Marquardt, grid search for rho_1 and rho, CRS
cesLmGrid2 <- cesEst( "yObs", xNames, data = MishraCES, rho1 = gridRhos[-1],
   rho = gridRhos[-c(12:13)], control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmGrid2 )
print( cesLmGrid2 )
summary( cesLmGrid2 )
coef( cesLmGrid2 )
vcov( cesLmGrid2 )
coef( summary( cesLmGrid2 ) )
fitted( cesLmGrid2 )
residuals( cesLmGrid2 )

## Levenberg-Marquardt, grid search for rho_1 and rho, VRS
cesLmVrsGrid2 <- cesEst( "yObs", xNames, data = MishraCES, vrs = TRUE,
   rho1 = gridRhos[-c(12:13)], rho = gridRhos[-1], 
   control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmVrsGrid2 )
print( cesLmVrsGrid2 )
summary( cesLmVrsGrid2 )
coef( cesLmVrsGrid2 )
vcov( cesLmVrsGrid2 )
coef( summary( cesLmVrsGrid2 ) )
fitted( cesLmVrsGrid2 )
residuals( cesLmVrsGrid2 )

## PORT, grid search for rho_1 and rho, CRS
cesPortGrid2 <- cesEst( "yObs", xNames, data = MishraCES, method = "PORT",
   rho1 = gridRhos[-c(12:13)], rho = gridRhos[-1], control = list( eval.max = 500 ) )
print.default( cesPortGrid2 )
print( cesPortGrid2 )
summary( cesPortGrid2 )
coef( cesPortGrid2 )
vcov( cesPortGrid2 )
coef( summary( cesPortGrid2 ) )
fitted( cesPortGrid2 )
residuals( cesPortGrid2 )

## PORT, grid search for rho_1 and rho, VRS
cesPortVrsGrid2 <- cesEst( "yObs", xNames, data = MishraCES, method = "PORT", 
   vrs = TRUE, rho1 = gridRhos[-1], rho = gridRhos[-c(12:13)], 
   control = list( eval.max = 500, iter.max = 500 ) )
print.default( cesPortVrsGrid2 )
print( cesPortVrsGrid2 )
summary( cesPortVrsGrid2 )
coef( cesPortVrsGrid2 )
vcov( cesPortVrsGrid2 )
coef( summary( cesPortVrsGrid2 ) )
fitted( cesPortVrsGrid2 )
residuals( cesPortVrsGrid2 )


