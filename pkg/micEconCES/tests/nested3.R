library( "micEconCES" )
library( "maxLik" )

data( "MishraCES" )

b <- c( "gamma_1" = 2, "gamma_2" = 200, "delta_1" = 0.6, "delta_2" = 0.7, 
   "rho_1" = 0.5, "rho" = 0.6 )
bVrs <- c( b, "nu" = 1.1 )
xNames <- c( "X1", "X2", "X3" )


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
yRho <- matrix( NA, nrow = length( rhos ), ncol = nrow( MishraCES ) )
rownames( yRho ) <- c( -(1:20), 0, (20:1) )

# calculate endogenous variables
coefRho <- b
for( i in 1:length( rhos ) ) {
   coefRho[ "rho" ] <- rhos[ i ]
   yRho[ i, ] <- cesCalc( xNames = xNames, data = MishraCES, coef = coefRho,
      nested = TRUE )
}

# print "raw" endogenous values
print( yRho )

# print endogenous variables for different rhos (adjusted with the y at rho=0)
for( i in 1:nrow( MishraCES ) ) {
   print( format( round( yRho[ , i, drop = FALSE ] - yRho[21,i], 11 ),
      scientific = FALSE ) )
}

## check cesCalc() with rho_1 equal to zero and close to zero
# matrix for returned endogenous variables
yRho1 <- matrix( NA, nrow = length( rhos ), ncol = nrow( MishraCES ) )
rownames( yRho1 ) <- c( -(1:20), 0, (20:1) )

# calculate endogenous variables
coefRho1 <- b
for( i in 1:length( rhos ) ) {
   coefRho1[ "rho_1" ] <- rhos[ i ]
   yRho1[ i, ] <- cesCalc( xNames = xNames, data = MishraCES, coef = coefRho1,
      nested = TRUE )
}

# print "raw" endogenous values
print( yRho1 )


## check cesCalc() with rho_1 and rho equal to zero and close to zero
# array for returned endogenous variables
yRho2 <- array( NA, c( length( rhos ), length( rhos ), 7 ) )
dimnames( yRho2 ) <- list( -20:20, -20:20, 1:7 )

# calculate endogenous variables
coefRho2 <- bVrs
for( i in 1:length( rhos ) ) {
   for( j in 1:length( rhos ) ) {
      coefRho2[ "rho" ] <- rhos[ i ]
      coefRho2[ "rho_1" ] <- rhos[ j ]
      yRho2[ i, j, ] <- cesCalc( xNames = xNames, data = MishraCES[ 1:7, ], 
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
# array for returned endogenous variables
deriv <- array( NA, c( length( rhos ), nrow( MishraCES ), length( bVrs ) ) )
dimnames( deriv ) <- list( rhos, 1:nrow( MishraCES ), names( bVrs ) )

coefRho <- bVrs
# calculate the derivatives
for( i in 1:length( rhos ) ) {
   # coefficients
   coefRho[ "rho" ] <- rhos[ i ]
   deriv[ i, , ] <- micEconCES:::cesDerivCoef( par = coefRho,
      xNames = xNames, data = MishraCES, nested = TRUE, vrs = TRUE,
      rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
}

# print array of derivatives
print( deriv )

# print derivatives for different rhos (adjusted with the derivatives at rho=0)
for( k in 1:dim( deriv )[3] ) {
   for( i in 1:nrow( MishraCES ) ) {
      print( format( round( deriv[ , i, k, drop = FALSE ] -
       deriv[ 21, i, k ], 11 ), scientific = FALSE ) )
   }
}

## check cesDerivCoef() with rho_1 equal to zero and close to zero
# array for returned endogenous variables
deriv1 <- array( NA, c( length( rhos ), nrow( MishraCES ), length( bVrs ) ) )
dimnames( deriv1 ) <- list( rhos, 1:nrow( MishraCES ), names( bVrs ) )

coefRho1 <- bVrs
# calculate the derivatives
for( i in 1:length( rhos ) ) {
   # coefficients
   coefRho1[ "rho_1" ] <- rhos[ i ]
   deriv1[ i, , ] <- micEconCES:::cesDerivCoef( par = coefRho1,
      xNames = xNames, data = MishraCES, nested = TRUE, vrs = TRUE,
      rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
}

# print array of derivatives
print( deriv1 )


## check cesDerivCoef() with rho and rho_1 equal to zero and close to zero
# array for returned endogenous variables
deriv2 <- array( NA, c( length( rhos ), length( rhos ), 7, length( bVrs ) ) )
dimnames( deriv2 ) <- list( c(-20:20), c(-20:20), 1:7, names( bVrs ) )

coefRho2 <- bVrs
# calculate the derivatives
for( i in 1:length( rhos ) ) {
   for( j in 1:length( rhos ) ) {
      # coefficients
      coefRho2[ "rho" ] <- rhos[ i ]
      coefRho2[ "rho_1" ] <- rhos[ j ]
      deriv2[ i, j, , ] <- micEconCES:::cesDerivCoef( par = coefRho2,
         xNames = xNames, data = MishraCES[1:7,], nested = TRUE, vrs = TRUE,
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


########## Grid Search for Rho ##############
## Levenberg-Marquardt, Grid Search, CRS
cesLmGrid <- cesEst( "yObs", xNames, data = MishraCES, rho = (1:12)/6-0.4,
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
   rho = (1:12)/6-0.4, control = nls.lm.control( maxiter = 200 ) )
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
   rho = (1:12)/6-0.4, control = list( eval.max = 500 ) )
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
   vrs = TRUE, rho = (1:12)/6-0.4, 
   control = list( eval.max = 500, iter.max = 500 ) )
print.default( cesPortVrsGrid )
print( cesPortVrsGrid )
summary( cesPortVrsGrid )
coef( cesPortVrsGrid )
vcov( cesPortVrsGrid )
coef( summary( cesPortVrsGrid ) )
fitted( cesPortVrsGrid )
residuals( cesPortVrsGrid )



