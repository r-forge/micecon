library( "micEconCES" )
library( maxLik )

options( max.print = 300000 )

data( "MishraCES" )

b <- c( "gamma" = 200 * 0.5^( 1 / 0.6 ), 
   "delta_1" = 0.6, "delta_2" = 0.3, "delta_3" = 0.5, "rho_1" = 0.5,
   "rho_2" = -0.17, "rho" = 0.6 )
bVrs <- c( b, "nu" = 1.05 )
bVrs[ "gamma" ] <- 200 * 0.5^( 1.05 / 0.6 )
xNames <- c( "X1", "X2", "X3", "X4" )


######################### checking cesCalc() ###############################
MishraCES$Y2 <- cesCalc( xNames = xNames, 
   data = MishraCES, coef = b, nested = TRUE )

MishraCES$Y2

all.equal( MishraCES$Y, MishraCES$Y2 )

# VRS
MishraCES$yVrs <- cesCalc( xNames = xNames, 
   data = MishraCES, coef = bVrs, nested = TRUE )
MishraCES$yVrs


## check cesCalc() with rho equal to zero and close to zero
# vector with values around 0 for coefficient "rho"
rhos <- c( -exp(-(1:20)),0,exp(-(20:1)) )

# matrix for returned endogenous variables
yRho <- matrix( NA, nrow = length( rhos ), ncol = nrow( MishraCES ) )
rownames( yRho ) <- c( -20:20 )

# calculate endogenous variables
coefRho <- bVrs
for( i in 1:length( rhos ) ) {
   coefRho[ "rho" ] <- rhos[ i ]
   yRho[ i, ] <- cesCalc( xNames = xNames, data = MishraCES, coef = coefRho,
      nested = TRUE )
}

# print "raw" endogenous values
print( yRho[ , 1:12 ] )


## check cesCalc() with rho_1 equal to zero and close to zero
# matrix for returned endogenous variables
yRho1 <- matrix( NA, nrow = length( rhos ), ncol = nrow( MishraCES ) )
rownames( yRho1 ) <- c( -20:20 )

# calculate endogenous variables
coefRho <- bVrs
for( i in 1:length( rhos ) ) {
   coefRho[ "rho_1" ] <- rhos[ i ]
   yRho1[ i, ] <- cesCalc( xNames = xNames, data = MishraCES, coef = coefRho,
      nested = TRUE )
}

# print "raw" endogenous values
print( yRho1[ , 1:12 ] )


## check cesCalc() with rho_2 equal to zero and close to zero
# matrix for returned endogenous variables
yRho2 <- matrix( NA, nrow = length( rhos ), ncol = nrow( MishraCES ) )
rownames( yRho2 ) <- c( -20:20 )

# calculate endogenous variables
coefRho <- bVrs
for( i in 1:length( rhos ) ) {
   coefRho[ "rho_2" ] <- rhos[ i ]
   yRho2[ i, ] <- cesCalc( xNames = xNames, data = MishraCES, coef = coefRho,
      nested = TRUE )
}

# print "raw" endogenous values
print( yRho2[ , 1:12 ] )


## check cesCalc() with rho_1 and rho equal to zero and close to zero
# array for returned endogenous variables
yRho10 <- array( NA, c( length( rhos ), length( rhos ), 7 ) )
dimnames( yRho10 ) <- list( -20:20, -20:20, 1:7 )

# calculate endogenous variables
coefRho <- bVrs
for( i in 1:length( rhos ) ) {
   for( j in 1:length( rhos ) ) {
      coefRho[ "rho_1" ] <- rhos[ i ]
      coefRho[ "rho" ] <- rhos[ j ]
      yRho10[ i, j, ] <- cesCalc( xNames = xNames, data = MishraCES[ 1:7, ], 
         coef = coefRho, nested = TRUE )
   }
}

# print "raw" endogenous values
print( yRho10 )


## check cesCalc() with rho_2 and rho equal to zero and close to zero
# array for returned endogenous variables
yRho20 <- array( NA, c( length( rhos ), length( rhos ), 7 ) )
dimnames( yRho20 ) <- list( -20:20, -20:20, 1:7 )

# calculate endogenous variables
coefRho <- bVrs
for( i in 1:length( rhos ) ) {
   for( j in 1:length( rhos ) ) {
      coefRho[ "rho_2" ] <- rhos[ i ]
      coefRho[ "rho" ] <- rhos[ j ]
      yRho20[ i, j, ] <- cesCalc( xNames = xNames, data = MishraCES[ 1:7, ], 
         coef = coefRho, nested = TRUE )
   }
}

# print "raw" endogenous values
print( yRho20 )


## check cesCalc() with rho_1 and rho_2 equal to zero and close to zero
# array for returned endogenous variables
yRho12 <- array( NA, c( length( rhos ), length( rhos ), 7 ) )
dimnames( yRho12 ) <- list( -20:20, -20:20, 1:7 )

# calculate endogenous variables
coefRho <- bVrs
for( i in 1:length( rhos ) ) {
   for( j in 1:length( rhos ) ) {
      coefRho[ "rho_1" ] <- rhos[ i ]
      coefRho[ "rho_2" ] <- rhos[ j ]
      yRho12[ i, j, ] <- cesCalc( xNames = xNames, data = MishraCES[ 1:7, ], 
         coef = coefRho, nested = TRUE )
   }
}

# print "raw" endogenous values
print( yRho12 )


## check cesCalc() with rho_1, rho_2, and rho equal to zero and close to zero
# array for returned endogenous variables
yRho120 <- array( NA, c( length( rhos ), length( rhos ), length( rhos ), 7 ) )
dimnames( yRho120 ) <- list( -20:20, -20:20, -20:20, 1:7 )

# calculate endogenous variables
coefRho <- bVrs
for( i in 1:length( rhos ) ) {
   for( j in 1:length( rhos ) ) {
      for( k in 1:length( rhos ) ) {
         coefRho[ "rho_1" ] <- rhos[ i ]
         coefRho[ "rho_2" ] <- rhos[ j ]
         coefRho[ "rho" ] <- rhos[ k ]
         yRho120[ i, j, k, ] <- cesCalc( xNames = xNames, 
            data = MishraCES[ 1:7, ], coef = coefRho, nested = TRUE )
      }
   }
}

# print "raw" endogenous values
print( yRho120[ , , c( 2, 16, 21, 26, 40 ), ] )
print( aperm( yRho120, c( 2, 3, 1, 4 ) )[ , , c( 2, 16, 21, 26, 40 ), ] )
print( aperm( yRho120, c( 3, 1, 2, 4 ) )[ , , c( 2, 16, 21, 26, 40 ), ] )


########################## checking cesDerivCoef ##############################
cesDeriv <- micEconCES:::cesDerivCoef( par = b, xNames = xNames, 
   data = MishraCES, vrs = FALSE, nested = TRUE, 
   rhoApprox = NULL )
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
   rhoApprox = NULL )
cesDerivVrsNum <- numericGradient( f, t0 = bVrs )
all.equal( cesDerivVrs, cesDerivVrsNum )
print( cesDerivVrs )


## check cesDerivCoef() with rho equal to zero and close to zero
# array for returned partial derivatives
deriv <- array( NA, c( length( rhos ), nrow( MishraCES ), length( bVrs ) ) )
dimnames( deriv ) <- list( -20:20, 1:nrow( MishraCES ), names( bVrs ) )

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
print( deriv[ , 1:12, ] )


## check cesDerivCoef() with rho_1 equal to zero and close to zero
# array for returned partial derivatives
deriv1 <- array( NA, c( length( rhos ), nrow( MishraCES ), length( bVrs ) ) )
dimnames( deriv1 ) <- list( -20:20, 1:nrow( MishraCES ), names( bVrs ) )

coefRho <- bVrs
# calculate the derivatives
for( i in 1:length( rhos ) ) {
   # coefficients
   coefRho[ "rho_1" ] <- rhos[ i ]
   deriv1[ i, , ] <- micEconCES:::cesDerivCoef( par = coefRho,
      xNames = xNames, data = MishraCES, nested = TRUE, vrs = TRUE,
      rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
}

# print array of derivatives
print( deriv1[ , 1:12, ] )


## check cesDerivCoef() with rho_2 equal to zero and close to zero
# array for returned partial derivatives
deriv2 <- array( NA, c( length( rhos ), nrow( MishraCES ), length( bVrs ) ) )
dimnames( deriv2 ) <- list( -20:20, 1:nrow( MishraCES ), names( bVrs ) )

coefRho <- bVrs
# calculate the derivatives
for( i in 1:length( rhos ) ) {
   # coefficients
   coefRho[ "rho_2" ] <- rhos[ i ]
   deriv2[ i, , ] <- micEconCES:::cesDerivCoef( par = coefRho,
      xNames = xNames, data = MishraCES, nested = TRUE, vrs = TRUE,
      rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
}

# print array of derivatives
print( deriv2[ , 1:12, ] )


## check cesDerivCoef() with rho_1 and rho equal to zero and close to zero
# array for returned partial derivatives
deriv10 <- array( NA, c( length( rhos ), length( rhos ), 3, length( bVrs ) ) )
dimnames( deriv10 ) <- list( c(-20:20), c(-20:20), 1:3, names( bVrs ) )

coefRho <- bVrs
# calculate the derivatives
for( i in 1:length( rhos ) ) {
   for( j in 1:length( rhos ) ) {
      # coefficients
      coefRho[ "rho_1" ] <- rhos[ i ]
      coefRho[ "rho" ] <- rhos[ j ]
      deriv10[ i, j, , ] <- micEconCES:::cesDerivCoef( par = coefRho,
         xNames = xNames, data = MishraCES[1:3,], nested = TRUE, vrs = TRUE,
         rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
   }
}

# print array of derivatives
print( deriv10 )


## check cesDerivCoef() with rho_2 and rho equal to zero and close to zero
# array for returned partial derivatives
deriv20 <- array( NA, c( length( rhos ), length( rhos ), 3, length( bVrs ) ) )
dimnames( deriv20 ) <- list( c(-20:20), c(-20:20), 1:3, names( bVrs ) )

coefRho <- bVrs
# calculate the derivatives
for( i in 1:length( rhos ) ) {
   for( j in 1:length( rhos ) ) {
      # coefficients
      coefRho[ "rho_2" ] <- rhos[ i ]
      coefRho[ "rho" ] <- rhos[ j ]
      deriv20[ i, j, , ] <- micEconCES:::cesDerivCoef( par = coefRho,
         xNames = xNames, data = MishraCES[1:3,], nested = TRUE, vrs = TRUE,
         rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
   }
}

# print array of derivatives
print( deriv20 )


## check cesDerivCoef() with rho_1 and rho_2 equal to zero and close to zero
# array for returned partial derivatives
deriv12 <- array( NA, c( length( rhos ), length( rhos ), 3, length( bVrs ) ) )
dimnames( deriv12 ) <- list( c(-20:20), c(-20:20), 1:3, names( bVrs ) )

coefRho <- bVrs
# calculate the derivatives
for( i in 1:length( rhos ) ) {
   for( j in 1:length( rhos ) ) {
      # coefficients
      coefRho[ "rho_1" ] <- rhos[ i ]
      coefRho[ "rho_2" ] <- rhos[ j ]
      deriv12[ i, j, , ] <- micEconCES:::cesDerivCoef( par = coefRho,
         xNames = xNames, data = MishraCES[1:3,], nested = TRUE, vrs = TRUE,
         rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
   }
}

# print array of derivatives
print( deriv12 )


## check cesDerivCoef() with rho_1, rho_2, and rho equal to zero and close to zero
# array for returned partial derivatives
deriv120 <- array( NA, 
   c( length( rhos ), length( rhos ), length( rhos ), 3, length( bVrs ) ) )
dimnames( deriv120 ) <- list( -20:20, -20:20, -20:20, 1:3, names( bVrs ) )

# calculate endogenous variables
coefRho <- bVrs
for( i in 1:length( rhos ) ) {
   for( j in 1:length( rhos ) ) {
      for( k in 1:length( rhos ) ) {
         coefRho[ "rho_1" ] <- rhos[ i ]
         coefRho[ "rho_2" ] <- rhos[ j ]
         coefRho[ "rho" ] <- rhos[ k ]
         deriv120[ i, j, k, , ] <- micEconCES:::cesDerivCoef( par = coefRho,
            xNames = xNames, data = MishraCES[1:3,], nested = TRUE, vrs = TRUE,
            rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
      }
   }
}

# print "raw" endogenous values
print( deriv120[ , , c( 2, 16, 21, 26, 40 ), , ] )
print( aperm( deriv120, c( 2, 3, 1, 4, 5 ) )[ , , c( 2, 16, 21, 26, 40 ), , ] )
print( aperm( deriv120, c( 3, 1, 2, 4, 5 ) )[ , , c( 2, 16, 21, 26, 40 ), , ] )


############################# checking cesEst #################################
set.seed( 345 )
MishraCES$yObs <- MishraCES$Y + 40 * rnorm( nrow( MishraCES ) )

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

## nls, CRS
cesNls <- cesEst( "yObs", xNames, data = MishraCES, method = "nls" )
print.default( cesNls )
print( cesNls )
summary( cesNls )
coef( cesNls )
vcov( cesNls )
coef( summary( cesNls ) )
fitted( cesNls )
residuals( cesNls )

## nls, VRS
cesNlsVrs <- cesEst( "yObs", xNames, data = MishraCES, method = "nls", 
   vrs = TRUE )
print.default( cesNlsVrs )
print( cesNlsVrs )
summary( cesNlsVrs )
coef( cesNlsVrs )
vcov( cesNlsVrs )
coef( summary( cesNlsVrs ) )
fitted( cesNlsVrs )
residuals( cesNlsVrs )


########## Estimation with Fixed Rhos ##############
## Levenberg-Marquardt, Fixed rho, CRS
cesLmR <- cesEst( "yObs", xNames, data = MishraCES,
   control = nls.lm.control( maxiter = 200 ), rho = 0.9 )
print.default( cesLmR )
print( cesLmR )
summary( cesLmR )
coef( cesLmR )
vcov( cesLmR )
coef( summary( cesLmR ) )
fitted( cesLmR )
residuals( cesLmR )

## Levenberg-Marquardt, Fixed rho, VRS
cesLmVrsR <- cesEst( "yObs", xNames, data = MishraCES,
   control = nls.lm.control( maxiter = 200 ), rho = 0 )
print.default( cesLmVrsR )
print( cesLmVrsR )
summary( cesLmVrsR )
coef( cesLmVrsR )
vcov( cesLmVrsR )
coef( summary( cesLmVrsR ) )
fitted( cesLmVrsR )
residuals( cesLmVrsR )

## Levenberg-Marquardt, Fixed rho_1, CRS
cesLmR1 <- cesEst( "yObs", xNames, data = MishraCES,
   control = nls.lm.control( maxiter = 200 ), rho1 = 0 )
print.default( cesLmR1 )
print( cesLmR1 )
summary( cesLmR1 )
coef( cesLmR1 )
vcov( cesLmR1 )
coef( summary( cesLmR1 ) )
fitted( cesLmR1 )
residuals( cesLmR1 )

## Levenberg-Marquardt, Fixed rho_1, VRS
cesLmVrsR1 <- cesEst( "yObs", xNames, data = MishraCES,
   control = nls.lm.control( maxiter = 200 ), rho1 = -0.1 )
print.default( cesLmVrsR1 )
print( cesLmVrsR1 )
summary( cesLmVrsR1 )
coef( cesLmVrsR1 )
vcov( cesLmVrsR1 )
coef( summary( cesLmVrsR1 ) )
fitted( cesLmVrsR1 )
residuals( cesLmVrsR1 )

## Levenberg-Marquardt, Fixed rho_2, CRS
cesLmR2 <- cesEst( "yObs", xNames, data = MishraCES,
   control = nls.lm.control( maxiter = 200 ), rho2 = 0.3 )
print.default( cesLmR2 )
print( cesLmR2 )
summary( cesLmR2 )
coef( cesLmR2 )
vcov( cesLmR2 )
coef( summary( cesLmR2 ) )
fitted( cesLmR2 )
residuals( cesLmR2 )

## Levenberg-Marquardt, Fixed rho_2, VRS
cesLmVrsR2 <- cesEst( "yObs", xNames, data = MishraCES,
   control = nls.lm.control( maxiter = 200 ), rho2 = 0 )
print.default( cesLmVrsR2 )
print( cesLmVrsR2 )
summary( cesLmVrsR2 )
coef( cesLmVrsR2 )
vcov( cesLmVrsR2 )
coef( summary( cesLmVrsR2 ) )
fitted( cesLmVrsR2 )
residuals( cesLmVrsR2 )

## Levenberg-Marquardt, Fixed rho1 and rho2, CRS
cesLmRR12 <- cesEst( "yObs", xNames, data = MishraCES,
   control = nls.lm.control( maxiter = 200 ), rho1 = 0.9, rho2 = 0 )
print.default( cesLmRR12 )
print( cesLmRR12 )
summary( cesLmRR12 )
coef( cesLmRR12 )
vcov( cesLmRR12 )
coef( summary( cesLmRR12 ) )
fitted( cesLmRR12 )
residuals( cesLmRR12 )

## Levenberg-Marquardt, Fixed rho and rho1, CRS
cesLmRR1 <- cesEst( "yObs", xNames, data = MishraCES,
   control = nls.lm.control( maxiter = 200 ), rho = 0.9, rho1 = 0 )
print.default( cesLmRR1 )
print( cesLmRR1 )
summary( cesLmRR1 )
coef( cesLmRR1 )
vcov( cesLmRR1 )
coef( summary( cesLmRR1 ) )
fitted( cesLmRR1 )
residuals( cesLmRR1 )

## Levenberg-Marquardt, Fixed rho and rho2, VRS
cesLmVrsRR2 <- cesEst( "yObs", xNames, data = MishraCES,
   control = nls.lm.control( maxiter = 200 ), rho = 0.2, rho2 = 0 )
print.default( cesLmVrsRR2 )
print( cesLmVrsRR2 )
summary( cesLmVrsRR2 )
coef( cesLmVrsRR2 )
vcov( cesLmVrsRR2 )
coef( summary( cesLmVrsRR2 ) )
fitted( cesLmVrsRR2 )
residuals( cesLmVrsRR2 )

## Levenberg-Marquardt, Fixed rho, rho1, and rho2, VRS
cesLmVrsRRR <- cesEst( "yObs", xNames, data = MishraCES,
   control = nls.lm.control( maxiter = 200 ), rho1 = 0.2, rho2 = 0.3, rho = -0.1 )
print.default( cesLmVrsRRR )
print( cesLmVrsRRR )
summary( cesLmVrsRRR )
coef( cesLmVrsRRR )
vcov( cesLmVrsRRR )
coef( summary( cesLmVrsRRR ) )
fitted( cesLmVrsRRR )
residuals( cesLmVrsRRR )


########## Grid Search for Rho_1, Rho_2, and Rho ##############
## Levenberg-Marquardt, Grid Search for rho, CRS
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

## Levenberg-Marquardt, Grid Search for rho, VRS
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

## PORT, Grid Search for rho, CRS
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

## PORT, Grid Search for rho, VRS
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

## Levenberg-Marquardt, Grid Search for rho_1, CRS
cesLmGrid1 <- cesEst( "yObs", xNames, data = MishraCES, rho1 = (1:8)/6,
   control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmGrid1 )
print( cesLmGrid1 )
summary( cesLmGrid1 )
coef( cesLmGrid1 )
vcov( cesLmGrid1 )
coef( summary( cesLmGrid1 ) )
fitted( cesLmGrid1 )
residuals( cesLmGrid1 )

## Levenberg-Marquardt, Grid Search for rho_2, VRS
cesLmVrsGrid2 <- cesEst( "yObs", xNames, data = MishraCES, vrs = TRUE,
   rho2 = (0:8)/6-0.6, control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmVrsGrid2 )
print( cesLmVrsGrid2 )
summary( cesLmVrsGrid2 )
coef( cesLmVrsGrid2 )
vcov( cesLmVrsGrid2 )
coef( summary( cesLmVrsGrid2 ) )
fitted( cesLmVrsGrid2 )
residuals( cesLmVrsGrid2 )

## Levenberg-Marquardt, Grid Search for rho_1 and rho_2, CRS
cesLmGrid12 <- cesEst( "yObs", xNames, data = MishraCES, rho1 = (1:8)/6,
   rho2 = 1:6/5-0.8, control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmGrid12 )
print( cesLmGrid12 )
summary( cesLmGrid12 )
coef( cesLmGrid12 )
vcov( cesLmGrid12 )
coef( summary( cesLmGrid12 ) )
fitted( cesLmGrid12 )
residuals( cesLmGrid12 )

## Levenberg-Marquardt, Grid Search for rho_1 and rho, VRS
cesLmVrsGrid10 <- cesEst( "yObs", xNames, data = MishraCES, vrs = TRUE,
   rho1 = (0:8)/6-0.1, rho = (1:7)/4-0.5,
   control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmVrsGrid10 )
print( cesLmVrsGrid10 )
summary( cesLmVrsGrid10 )
coef( cesLmVrsGrid10 )
vcov( cesLmVrsGrid10 )
coef( summary( cesLmVrsGrid10 ) )
fitted( cesLmVrsGrid10 )
residuals( cesLmVrsGrid10 )

## Levenberg-Marquardt, Grid Search for rho_2 and rho, CRS
cesLmGrid20 <- cesEst( "yObs", xNames, data = MishraCES, rho2 = (1:8)/7-0.7,
   rho = 1:6/7-0.1, control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmGrid20 )
print( cesLmGrid20 )
summary( cesLmGrid20 )
coef( cesLmGrid20 )
vcov( cesLmGrid20 )
coef( summary( cesLmGrid20 ) )
fitted( cesLmGrid20 )
residuals( cesLmGrid20 )

## Levenberg-Marquardt, Grid Search for rho_1, rho_2, and rho, VRS
cesLmVrsGrid123 <- cesEst( "yObs", xNames, data = MishraCES, vrs = TRUE,
   rho1 = (0:5)/5+0.3, rho2 = (0:4)/3-0.9, rho = (0:5)/5+0.1,
   control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmVrsGrid123 )
print( cesLmVrsGrid123 )
summary( cesLmVrsGrid123 )
coef( cesLmVrsGrid123 )
vcov( cesLmVrsGrid123 )
coef( summary( cesLmVrsGrid123 ) )
fitted( cesLmVrsGrid123 )
residuals( cesLmVrsGrid123 )


