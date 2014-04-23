library( "micEconCES" )
library( maxLik )

options( max.print = 300000 )

data( "MishraCES" )
MishraCES$time <- c( 0:49 )

b <- c( "gamma" = 200 * 0.5^( 1 / 0.6 ), 
   "delta_1" = 0.6, "delta_2" = 0.3, "delta" = 0.5, "rho_1" = 0.5,
   "rho_2" = -0.17, "rho" = 0.6 )
bVrs <- c( b, "nu" = 1.05 )
bVrs[ "gamma" ] <- 200 * 0.5^( 1.05 / 0.6 )
xNames <- c( "X1", "X2", "X3", "X4" )

# with technological change
bTc <- c( b[ 1 ], lambda = 0.02, b[ -1 ] )
bTcVrs <- c( bVrs[ 1 ], lambda = 0.02, bVrs[ -1 ] )


######################### checking cesCalc() ###############################
MishraCES$Y2 <- cesCalc( xNames = xNames, 
   data = MishraCES, coef = b, nested = TRUE )

MishraCES$Y2

all.equal( MishraCES$Y, MishraCES$Y2 )

# VRS
MishraCES$yVrs <- cesCalc( xNames = xNames, 
   data = MishraCES, coef = bVrs, nested = TRUE )
MishraCES$yVrs

# with technological change
MishraCES$yTc <- cesCalc( xNames = xNames, tName = "time",
   data = MishraCES, coef = bTc, nested = TRUE )
MishraCES$yTc
all.equal( MishraCES$yTc, MishraCES$Y2 * exp( bTc[ "lambda" ] * 0:49 ) )
# check if removing the names of the coefficients makes a difference
all.equal( MishraCES$yTc,
   cesCalc( xNames = xNames, tName = "time", data = MishraCES, 
      coef = unname( bTc ), nested = TRUE ) )
# check if permuting the coefficients makes a difference
all.equal( MishraCES$yTc,
   cesCalc( xNames = xNames, tName = "time", data = MishraCES, 
      coef = sample( bTc, 8 ), nested = TRUE ) )

# with technological change and VRS
MishraCES$yTcVrs <- cesCalc( xNames = xNames, tName = "time", 
   data = MishraCES, coef = bTcVrs, nested = TRUE )
MishraCES$yTcVrs
all.equal( MishraCES$yTcVrs, MishraCES$yVrs * exp( bTcVrs[ "lambda" ] * 0:49 ) )
# check if removing the names of the coefficients makes a difference
all.equal( MishraCES$yTcVrs,
   cesCalc( xNames = xNames, tName = "time", data = MishraCES, 
      coef = unname( bTcVrs ), nested = TRUE ) )
# check if permuting the coefficients makes a difference
all.equal( MishraCES$yTcVrs,
   cesCalc( xNames = xNames, tName = "time", data = MishraCES, 
      coef = sample( bTcVrs, 9 ), nested = TRUE ) )


## check cesCalc() with rho equal to zero and close to zero
# vector with values around 0 for coefficient "rho"
rhos <- c( -exp(-(1:20)),0,exp(-(20:1)) )

# matrix for returned endogenous variables
yRho <- matrix( NA, nrow = length( rhos ), ncol = 7 )
rownames( yRho ) <- c( -20:20 )

# calculate endogenous variables
coefRho <- bTcVrs
for( i in 1:length( rhos ) ) {
   coefRho[ "rho" ] <- rhos[ i ]
   yRho[ i, ] <- cesCalc( xNames = xNames, tName = "time",
      data = MishraCES[1:7,], coef = coefRho,
      nested = TRUE )
}

# print "raw" endogenous values
print( yRho )


## check cesCalc() with rho_1 equal to zero and close to zero
# matrix for returned endogenous variables
yRho1 <- matrix( NA, nrow = length( rhos ), ncol = 7 )
rownames( yRho1 ) <- c( -20:20 )

# calculate endogenous variables
coefRho <- bTcVrs
for( i in 1:length( rhos ) ) {
   coefRho[ "rho_1" ] <- rhos[ i ]
   yRho1[ i, ] <- cesCalc( xNames = xNames, tName = "time",
      data = MishraCES[1:7, ], coef = coefRho,
      nested = TRUE )
}

# print "raw" endogenous values
print( yRho1 )


## check cesCalc() with rho_2 equal to zero and close to zero
# matrix for returned endogenous variables
yRho2 <- matrix( NA, nrow = length( rhos ), ncol = 7 )
rownames( yRho2 ) <- c( -20:20 )

# calculate endogenous variables
coefRho <- bTcVrs
for( i in 1:length( rhos ) ) {
   coefRho[ "rho_2" ] <- rhos[ i ]
   yRho2[ i, ] <- cesCalc( xNames = xNames, tName = "time", 
      data = MishraCES[1:7,], coef = coefRho,
      nested = TRUE )
}

# print "raw" endogenous values
print( yRho2 )


## check cesCalc() with rho_1 and rho equal to zero and close to zero
# array for returned endogenous variables
yRho10 <- array( NA, c( length( rhos ), length( rhos ), 2 ) )
dimnames( yRho10 ) <- list( -20:20, -20:20, 1:2 )

# calculate endogenous variables
coefRho <- bTcVrs
for( i in 1:length( rhos ) ) {
   for( j in 1:length( rhos ) ) {
      coefRho[ "rho_1" ] <- rhos[ i ]
      coefRho[ "rho" ] <- rhos[ j ]
      yRho10[ i, j, ] <- cesCalc( xNames = xNames, tName = "time",
         data = MishraCES[ 1:2, ], 
         coef = coefRho, nested = TRUE )
   }
}

# print "raw" endogenous values
print( yRho10 )


## check cesCalc() with rho_2 and rho equal to zero and close to zero
# array for returned endogenous variables
yRho20 <- array( NA, c( length( rhos ), length( rhos ), 2 ) )
dimnames( yRho20 ) <- list( -20:20, -20:20, 1:2 )

# calculate endogenous variables
coefRho <- bTcVrs
for( i in 1:length( rhos ) ) {
   for( j in 1:length( rhos ) ) {
      coefRho[ "rho_2" ] <- rhos[ i ]
      coefRho[ "rho" ] <- rhos[ j ]
      yRho20[ i, j, ] <- cesCalc( xNames = xNames, tName = "time",
         data = MishraCES[ 1:2, ], 
         coef = coefRho, nested = TRUE )
   }
}

# print "raw" endogenous values
print( yRho20 )


## check cesCalc() with rho_1 and rho_2 equal to zero and close to zero
# array for returned endogenous variables
yRho12 <- array( NA, c( length( rhos ), length( rhos ), 2 ) )
dimnames( yRho12 ) <- list( -20:20, -20:20, 1:2 )

# calculate endogenous variables
coefRho <- bTcVrs
for( i in 1:length( rhos ) ) {
   for( j in 1:length( rhos ) ) {
      coefRho[ "rho_1" ] <- rhos[ i ]
      coefRho[ "rho_2" ] <- rhos[ j ]
      yRho12[ i, j, ] <- cesCalc( xNames = xNames, tName = "time",
         data = MishraCES[ 1:2, ], 
         coef = coefRho, nested = TRUE )
   }
}

# print "raw" endogenous values
print( yRho12 )


## check cesCalc() with rho_1, rho_2, and rho equal to zero and close to zero
# array for returned endogenous variables
rhos2 <- rhos[ seq( 1, length( rhos ), 2 ) ]
rhos2Num <- (-20:20)[ seq( 1, length( rhos ), 2 ) ]
yRho120 <- array( NA, c( length( rhos2 ), length( rhos2 ), length( rhos2 ), 2 ) )
dimnames( yRho120 ) <- list( rhos2Num, rhos2Num, rhos2Num, 1:2 )

# calculate endogenous variables
coefRho <- bTcVrs
for( i in 1:length( rhos2 ) ) {
   for( j in 1:length( rhos2 ) ) {
      for( k in 1:length( rhos2 ) ) {
         coefRho[ "rho_1" ] <- rhos2[ i ]
         coefRho[ "rho_2" ] <- rhos2[ j ]
         coefRho[ "rho" ] <- rhos2[ k ]
         yRho120[ i, j, k, ] <- cesCalc( xNames = xNames, tName = "time",
            data = MishraCES[ 1:2, ], coef = coefRho, nested = TRUE )
      }
   }
}

# print "raw" endogenous values
print( yRho120[ , , c( 2, 8, 11, 14, 20 ), ] )
print( aperm( yRho120, c( 2, 3, 1, 4 ) )[ , , c( 2, 8, 11, 14, 20 ), ] )
print( aperm( yRho120, c( 3, 1, 2, 4 ) )[ , , c( 2, 8, 11, 14, 20 ), ] )


########################## checking cesDerivCoef ##############################
cesDeriv <- micEconCES:::cesDerivCoef( par = bTc, xNames = xNames, 
   tName = "time", data = MishraCES, vrs = FALSE, nested = TRUE, 
   rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
f <- function( par ) {
   return( cesCalc( xNames = xNames, tName = "time", data = MishraCES, 
      coef = par, nested = TRUE ) )
}
cesDerivNum <- numericGradient( f, t0 = bTc )
all.equal( cesDeriv, cesDerivNum )
print( cesDeriv )

# VRS
cesDerivVrs <- micEconCES:::cesDerivCoef( par = bTcVrs, xNames = xNames, 
   tName = "time", data = MishraCES, vrs = TRUE, nested = TRUE,
   rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
cesDerivVrsNum <- numericGradient( f, t0 = bTcVrs )
all.equal( cesDerivVrs, cesDerivVrsNum )
print( cesDerivVrs )


## check cesDerivCoef() with rho equal to zero and close to zero
# array for returned partial derivatives
deriv <- array( NA, c( length( rhos ), 7, length( bTcVrs ) ) )
dimnames( deriv ) <- list( -20:20, 1:7, names( bTcVrs ) )

coefRho <- bTcVrs
# calculate the derivatives
for( i in 1:length( rhos ) ) {
   # coefficients
   coefRho[ "rho" ] <- rhos[ i ]
   deriv[ i, , ] <- micEconCES:::cesDerivCoef( par = coefRho, xNames = xNames, 
      tName = "time", data = MishraCES[1:7,], nested = TRUE, vrs = TRUE,
      rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
}

# print array of derivatives
print( deriv )


## check cesDerivCoef() with rho_1 equal to zero and close to zero
# array for returned partial derivatives
deriv1 <- array( NA, c( length( rhos ), 7, length( bTcVrs ) ) )
dimnames( deriv1 ) <- list( -20:20, 1:7, names( bTcVrs ) )

coefRho <- bTcVrs
# calculate the derivatives
for( i in 1:length( rhos ) ) {
   # coefficients
   coefRho[ "rho_1" ] <- rhos[ i ]
   deriv1[ i, , ] <- micEconCES:::cesDerivCoef( par = coefRho, xNames = xNames, 
      tName = "time", data = MishraCES[1:7,], nested = TRUE, vrs = TRUE,
      rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
}

# print array of derivatives
print( deriv1 )


## check cesDerivCoef() with rho_2 equal to zero and close to zero
# array for returned partial derivatives
deriv2 <- array( NA, c( length( rhos ), 7, length( bTcVrs ) ) )
dimnames( deriv2 ) <- list( -20:20, 1:7, names( bTcVrs ) )

coefRho <- bTcVrs
# calculate the derivatives
for( i in 1:length( rhos ) ) {
   # coefficients
   coefRho[ "rho_2" ] <- rhos[ i ]
   deriv2[ i, , ] <- micEconCES:::cesDerivCoef( par = coefRho, xNames = xNames, 
      tName = "time", data = MishraCES[1:7,], nested = TRUE, vrs = TRUE,
      rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
}

# print array of derivatives
print( deriv2 )


## check cesDerivCoef() with rho_1 and rho equal to zero and close to zero
# array for returned partial derivatives
deriv10 <- array( NA, c( length( rhos ), length( rhos ), 2, length( bTcVrs ) ) )
dimnames( deriv10 ) <- list( c(-20:20), c(-20:20), 1:2, names( bTcVrs ) )

coefRho <- bTcVrs
# calculate the derivatives
for( i in 1:length( rhos ) ) {
   for( j in 1:length( rhos ) ) {
      # coefficients
      coefRho[ "rho_1" ] <- rhos[ i ]
      coefRho[ "rho" ] <- rhos[ j ]
      deriv10[ i, j, , ] <- micEconCES:::cesDerivCoef( par = coefRho,
         xNames = xNames, tName = "time", data = MishraCES[1:2,], 
         nested = TRUE, vrs = TRUE,
         rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
   }
}

# print array of derivatives
print( deriv10 )


## check cesDerivCoef() with rho_2 and rho equal to zero and close to zero
# array for returned partial derivatives
deriv20 <- array( NA, c( length( rhos ), length( rhos ), 2, length( bTcVrs ) ) )
dimnames( deriv20 ) <- list( c(-20:20), c(-20:20), 1:2, names( bTcVrs ) )

coefRho <- bTcVrs
# calculate the derivatives
for( i in 1:length( rhos ) ) {
   for( j in 1:length( rhos ) ) {
      # coefficients
      coefRho[ "rho_2" ] <- rhos[ i ]
      coefRho[ "rho" ] <- rhos[ j ]
      deriv20[ i, j, , ] <- micEconCES:::cesDerivCoef( par = coefRho,
         xNames = xNames, tName = "time", data = MishraCES[1:2,], 
         nested = TRUE, vrs = TRUE,
         rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
   }
}

# print array of derivatives
print( deriv20 )


## check cesDerivCoef() with rho_1 and rho_2 equal to zero and close to zero
# array for returned partial derivatives
deriv12 <- array( NA, c( length( rhos ), length( rhos ), 2, length( bTcVrs ) ) )
dimnames( deriv12 ) <- list( c(-20:20), c(-20:20), 1:2, names( bTcVrs ) )

coefRho <- bTcVrs
# calculate the derivatives
for( i in 1:length( rhos ) ) {
   for( j in 1:length( rhos ) ) {
      # coefficients
      coefRho[ "rho_1" ] <- rhos[ i ]
      coefRho[ "rho_2" ] <- rhos[ j ]
      deriv12[ i, j, , ] <- micEconCES:::cesDerivCoef( par = coefRho,
         xNames = xNames, tName = "time", data = MishraCES[1:2,], 
         nested = TRUE, vrs = TRUE,
         rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
   }
}

# print array of derivatives
print( deriv12 )


## check cesDerivCoef() with rho_1, rho_2, and rho equal to zero and close to zero
# array for returned partial derivatives
rhos3 <- rhos[ seq( 3, length( rhos ), 3 ) ]
rhos3Nums <- (-20:20)[ seq( 3, length( rhos ), 3 ) ]
deriv120 <- array( NA, 
   c( length( rhos3 ), length( rhos3 ), length( rhos3 ), 2, length( bTcVrs ) ) )
dimnames( deriv120 ) <- list( rhos3Nums, rhos3Nums, rhos3Nums, 1:2, names( bTcVrs ) )

# calculate endogenous variables
coefRho <- bTcVrs
for( i in 1:length( rhos3 ) ) {
   for( j in 1:length( rhos3 ) ) {
      for( k in 1:length( rhos3 ) ) {
         coefRho[ "rho_1" ] <- rhos3[ i ]
         coefRho[ "rho_2" ] <- rhos3[ j ]
         coefRho[ "rho" ] <- rhos3[ k ]
         deriv120[ i, j, k, , ] <- micEconCES:::cesDerivCoef( par = coefRho,
            xNames = xNames, tName = "time", data = MishraCES[1:2,], 
            nested = TRUE, vrs = TRUE,
            rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
      }
   }
}

# print "raw" endogenous values
print( deriv120[ , , c( 1, 5, 7, 9, 13 ), , ] )
print( aperm( deriv120, c( 2, 3, 1, 4, 5 ) )[ , , c( 1, 5, 7, 9, 13 ), , ] )
print( aperm( deriv120, c( 3, 1, 2, 4, 5 ) )[ , , c( 1, 5, 7, 9, 13 ), , ] )


############################# checking cesEst #################################
set.seed( 345 )
options( digits = 3 )
MishraCES$yObs <- MishraCES$Y + 40 * rnorm( nrow( MishraCES ) )
MishraCES$yTcObs <- MishraCES$yTc + 40 * rnorm( nrow( MishraCES ) )
MishraCES$yTcVrsObs <- MishraCES$yTcVrs + 40 * rnorm( nrow( MishraCES ) )
MishraCES$yMeObs <- MishraCES$Y * exp( 0.3 * rnorm( nrow( MishraCES ) ) )
MishraCES$yTcMeObs <- MishraCES$yTc * exp( 0.3 * rnorm( nrow( MishraCES ) ) )
MishraCES$yTcMeVrsObs <- MishraCES$yTcVrs * exp( 0.3 * rnorm( nrow( MishraCES ) ) )

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

## Nelder-Mead, TC, CRS
cesNmTc <- cesEst( "yTcObs", xNames, tName = "time", data = MishraCES, 
   method = "Nelder-Mead", control = list( maxit = 2000 )  )
print.default( cesNmTc ) 
print( cesNmTc )
summary( cesNmTc )

## Nelder-Mead, TC, VRS
cesNmTcVrs <- cesEst( "yTcVrsObs", xNames, tName = "time", data = MishraCES, 
   vrs = TRUE, method = "NM", control = list( maxit = 5000 ) )
print.default( cesNmTcVrs )
print( cesNmTcVrs )
summary( cesNmTcVrs )

## Nelder-Mead, TC, multErr, CRS
cesNmTcMe <- cesEst( "yTcMeObs", xNames, tName = "time", data = MishraCES, 
   method = "Nelder-Mead", multErr = TRUE, control = list( maxit = 2000 ) )
print.default( cesNmTcMe ) 
print( cesNmTcMe )
summary( cesNmTcMe )

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

## Conjugate Gradients, TC, CRS
cesCgTc <- cesEst( "yTcVrsObs", xNames, tName = "time", data = MishraCES, 
   method = "CG" )
print.default( cesCgTc )
print( cesCgTc )
summary( cesCgTc )

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

## Simulated Annealing, TC, CRS
cesSannTc <- cesEst( "yTcObs", xNames, tName = "time", data = MishraCES, 
   method = "SANN" )
print.default( cesSannTc )
print( cesSannTc )
summary( cesSannTc )

## Simulated Annealing, multErr, VRS
cesSannMeVrs <- cesEst( "yMeObs", xNames, data = MishraCES, method = "SANN", 
   vrs = TRUE, multErr = TRUE )
print.default( cesSannMeVrs )
print( cesSannMeVrs )
summary( cesSannMeVrs )
vcov( cesSannMeVrs )

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

## BFGS, TC, CRS
cesBfgsTc <- cesEst( "yTcObs", xNames, tName = "time", data = MishraCES, 
   method = "BFGS", control = list( maxit = 500 ) )
print.default( cesBfgsTc )
print( cesBfgsTc )
summary( cesBfgsTc )

## BFGS, TC, VRS
cesBfgsTcVrs <- cesEst( "yTcVrsObs", xNames, tName = "time", data = MishraCES, 
   method = "BFGS", vrs = TRUE, control = list( maxit = 500 ) )
print.default( cesBfgsTcVrs )
print( cesBfgsTcVrs )
summary( cesBfgsTcVrs )

## BFGS, CRS, multErr
cesBfgsMe <- cesEst( "yMeObs", xNames, data = MishraCES, method = "BFGS",
   multErr = TRUE, control = list( maxit = 500 ) )
print.default( cesBfgsMe )
print( cesBfgsMe )
summary( cesBfgsMe )
vcov( cesBfgsMe )

## L-BFGS-B with constrained parameters, CRS
cesBfgsCon <- cesEst( "yObs", xNames, data = MishraCES, method = "L-BFGS-B" )
print.default( cesBfgsCon )
print( cesBfgsCon )
summary( cesBfgsCon )
summary( cesBfgsCon, ela = FALSE )
print( summary( cesBfgsCon ), ela = FALSE )
summary( cesBfgsCon )$elaCov
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

## L-BFGS-B with constrained parameters, TC, CRS
cesBfgsConTc <- cesEst( "yTcObs", xNames, tName = "time", data = MishraCES, 
   method = "L-BFGS-B", control = list( maxit = 500 ) )
print.default( cesBfgsConTc )
print( cesBfgsConTc )
summary( cesBfgsConTc )

## L-BFGS-B with constrained parameters, TC, VRS
cesBfgsConTcVrs <- cesEst( "yTcVrsObs", xNames, tName = "time", data = MishraCES, 
   method = "L-BFGS-B", vrs = TRUE, control = list( maxit = 500 ) )
print.default( cesBfgsConTcVrs )
print( cesBfgsConTcVrs )
summary( cesBfgsConTcVrs )

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

## Levenberg-Marquardt, TC, CRS
cesLmTc <- cesEst( "yTcObs", xNames, tName = "time", data = MishraCES,
   control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmTc )
print( cesLmTc )
summary( cesLmTc )

## Levenberg-Marquardt, TC, VRS
cesLmTcVrs <- cesEst( "yTcVrsObs", xNames, tName = "time", data = MishraCES, 
   vrs = TRUE, control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmTcVrs )
print( cesLmTcVrs )
summary( cesLmTcVrs )

## Levenberg-Marquardt, TC, multErr, VRS
cesLmTcMeVrs <- cesEst( "yTcMeVrsObs", xNames, tName = "time", data = MishraCES, 
   vrs = TRUE, multErr = TRUE, control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmTcMeVrs )
print( cesLmTcMeVrs )
summary( cesLmTcMeVrs )

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

## Newton-type, TC, CRS
cesNewtonTc <- cesEst( "yTcObs", xNames, tName = "time", data = MishraCES, 
   method = "Newton", iterlim = 500 )
print.default( cesNewtonTc )
print( cesNewtonTc )
summary( cesNewtonTc )

## Newton-type, TC, VRS
cesNewtonTcVrs <- cesEst( "yTcVrsObs", xNames, tName = "time", data = MishraCES, 
   method = "Newton", vrs = TRUE, iterlim = 500 )
print.default( cesNewtonTcVrs )
print( cesNewtonTcVrs )
summary( cesNewtonTcVrs )

## Newton-type, TC, multErr, VRS
cesNewtonTcMeVrs <- cesEst( "yTcMeVrsObs", xNames, tName = "time", data = MishraCES, 
   method = "Newton", vrs = TRUE, multErr = TRUE, iterlim = 200 )
print.default( cesNewtonTcMeVrs )
print( cesNewtonTcMeVrs )
summary( cesNewtonTcMeVrs )

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

## PORT, TC, CRS
cesPortTc <- cesEst( "yTcObs", xNames, tName = "time", data = MishraCES, 
   method = "PORT", control = list( eval.max = 500 ) )
print.default( cesPortTc )
print( cesPortTc )
summary( cesPortTc )

## PORT, TC, VRS
cesPortTcVrs <- cesEst( "yTcVrsObs", xNames, tName = "time", data = MishraCES, 
   method = "PORT", vrs = TRUE, 
   control = list( eval.max = 500, iter.max = 500 ) )
print.default( cesPortTcVrs )
print( cesPortTcVrs )
summary( cesPortTcVrs )

## PORT, TC, multErr, VRS
cesPortTcMeVrs <- cesEst( "yTcMeVrsObs", xNames, tName = "time", data = MishraCES, 
   method = "PORT", vrs = TRUE, multErr = TRUE,
   control = list( eval.max = 500, iter.max = 500 ) )
print.default( cesPortTcMeVrs )
print( cesPortTcMeVrs )
summary( cesPortTcMeVrs )

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

## DE, TC, CRS
cesDeTc <- cesEst( "yTcObs", xNames, tName = "time", data = MishraCES, 
   method = "DE", control = DEoptim.control( trace = FALSE, NP = 80 ) )
print.default( cesDeTc )
print( cesDeTc )
summary( cesDeTc )

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

## nls, TC, VRS
cesNlsTcVrs <- cesEst( "yTcVrsObs", xNames, tName = "time", data = MishraCES, 
   method = "nls", vrs = TRUE )
print.default( cesNlsTcVrs )
print( cesNlsTcVrs )
summary( cesNlsTcVrs )

## nls, TC, multErr, VRS
cesNlsTcMeVrs <- cesEst( "yTcMeVrsObs", xNames, tName = "time", data = MishraCES, 
   method = "nls", vrs = TRUE, multErr = TRUE )
print.default( cesNlsTcMeVrs )
print( cesNlsTcMeVrs )
summary( cesNlsTcMeVrs )


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

## Levenberg-Marquardt, Fixed rho, TC, CRS
cesLmTcR <- cesEst( "yTcObs", xNames, tName = "time", data = MishraCES,
   control = nls.lm.control( maxiter = 200 ), rho = 0.9 )
print.default( cesLmTcR )
print( cesLmTcR )
summary( cesLmTcR )

## Levenberg-Marquardt, Fixed rho at 0, VRS
cesLmVrsR <- cesEst( "yObs", xNames, data = MishraCES, vrs = TRUE,
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
cesLmVrsR1 <- cesEst( "yObs", xNames, data = MishraCES, vrs = TRUE,
   control = nls.lm.control( maxiter = 200 ), rho1 = -0.1 )
print.default( cesLmVrsR1 )
print( cesLmVrsR1 )
summary( cesLmVrsR1 )
coef( cesLmVrsR1 )
vcov( cesLmVrsR1 )
coef( summary( cesLmVrsR1 ) )
fitted( cesLmVrsR1 )
residuals( cesLmVrsR1 )

## Levenberg-Marquardt, Fixed rho_1, TC, VRS
cesLmTcVrsR1 <- cesEst( "yTcVrsObs", xNames, tName = "time", data = MishraCES,
   vrs = TRUE, control = nls.lm.control( maxiter = 200 ), rho1 = -0.1 )
print.default( cesLmTcVrsR1 )
print( cesLmTcVrsR1 )
summary( cesLmTcVrsR1 )

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

## Levenberg-Marquardt, Fixed rho_2, TC, CRS
cesLmTcR2 <- cesEst( "yTcObs", xNames, tName = "time", data = MishraCES,
   control = nls.lm.control( maxiter = 200 ), rho2 = 0.3 )
print.default( cesLmTcR2 )
print( cesLmTcR2 )
summary( cesLmTcR2 )

## Levenberg-Marquardt, Fixed rho_2, VRS
cesLmVrsR2 <- cesEst( "yObs", xNames, data = MishraCES, vrs = TRUE,
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
cesLmVrsRR2 <- cesEst( "yObs", xNames, data = MishraCES, vrs = TRUE,
   control = nls.lm.control( maxiter = 200 ), rho = 0.2, rho2 = 0 )
print.default( cesLmVrsRR2 )
print( cesLmVrsRR2 )
summary( cesLmVrsRR2 )
coef( cesLmVrsRR2 )
vcov( cesLmVrsRR2 )
coef( summary( cesLmVrsRR2 ) )
fitted( cesLmVrsRR2 )
residuals( cesLmVrsRR2 )

## Levenberg-Marquardt, Fixed rho and rho2, TC, VRS
cesLmTcVrsRR2 <- cesEst( "yTcVrsObs", xNames, tName = "time", data = MishraCES,
   vrs = TRUE, control = nls.lm.control( maxiter = 200 ), rho = 0.2, rho2 = 0 )
print.default( cesLmTcVrsRR2 )
print( cesLmTcVrsRR2 )
summary( cesLmTcVrsRR2 )

## Levenberg-Marquardt, Fixed rho, rho1, and rho2, VRS
cesLmVrsRRR <- cesEst( "yObs", xNames, data = MishraCES, vrs = TRUE,
   control = nls.lm.control( maxiter = 200 ), rho1 = 0.2, rho2 = 0.3, rho = -0.1 )
print.default( cesLmVrsRRR )
print( cesLmVrsRRR )
summary( cesLmVrsRRR )
coef( cesLmVrsRRR )
vcov( cesLmVrsRRR )
coef( summary( cesLmVrsRRR ) )
fitted( cesLmVrsRRR )
residuals( cesLmVrsRRR )

## Levenberg-Marquardt, Fixed rho, rho1, and rho2, TC, VRS
cesLmTcVrsRRR <- cesEst( "yTcVrsObs", xNames, tName = "time", data = MishraCES,
   vrs = TRUE, control = nls.lm.control( maxiter = 200 ), 
   rho1 = 0.2, rho2 = 0.3, rho = -0.1 )
print.default( cesLmTcVrsRRR )
print( cesLmTcVrsRRR )
summary( cesLmTcVrsRRR )

## Levenberg-Marquardt, Fixed rho, rho1, and rho2, multErr, VRS
cesLmMeVrsRRR <- cesEst( "yMeObs", xNames, data = MishraCES, vrs = TRUE,
   multErr = TRUE, control = nls.lm.control( maxiter = 200 ), 
   rho1 = 0.2, rho2 = 0.3, rho = -0.1 )
print.default( cesLmMeVrsRRR )
print( cesLmMeVrsRRR )
summary( cesLmMeVrsRRR )
vcov( cesLmMeVrsRRR )


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

## Levenberg-Marquardt, Grid Search for rho, TC, CRS
cesLmGridTc <- cesEst( "yObs", xNames, tName = "time", data = MishraCES, 
   rho = (1:6)/3-0.5, control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmGridTc )
print( cesLmGridTc )
summary( cesLmGridTc )

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

## PORT, Grid Search for rho, TC, VRS
cesPortTcVrsGrid1 <- cesEst( "yObs", xNames, tName = "time", data = MishraCES, 
   method = "PORT", vrs = TRUE, rho1 = (1:6)/3-0.5, 
   control = list( eval.max = 500, iter.max = 500 ) )
print.default( cesPortTcVrsGrid1 )
print( cesPortTcVrsGrid1 )
summary( cesPortTcVrsGrid1 )

## PORT, Grid Search for rho, TC, multErr CRS
cesPortTcMeGrid <- cesEst( "yTcMeObs", xNames, tName = "time", data = MishraCES, 
   method = "PORT", multErr = TRUE, rho = (1:12)/6-0.4, 
   control = list( eval.max = 500 ) )
print.default( cesPortTcMeGrid )
print( cesPortTcMeGrid )
summary( cesPortTcMeGrid )
vcov( cesPortTcMeGrid )

## Levenberg-Marquardt, Grid Search for rho_2, VRS
cesLmVrsGrid2 <- cesEst( "yObs", xNames, data = MishraCES, vrs = TRUE,
   rho2 = (0:8)/6-0.6, control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmVrsGrid2 )
print( cesLmVrsGrid2 )
summary( cesLmVrsGrid2 )
summary( cesLmVrsGrid2, ela = FALSE )
print( summary( cesLmVrsGrid2 ), ela = FALSE )
summary( cesLmVrsGrid2 )$elaCov
coef( cesLmVrsGrid2 )
vcov( cesLmVrsGrid2 )
coef( summary( cesLmVrsGrid2 ) )
fitted( cesLmVrsGrid2 )
residuals( cesLmVrsGrid2 )

## BFGS, Grid Search for rho_2, TC, VRS
cesBfgsVrsGrid2 <- cesEst( "yObs", xNames, tName = "time", data = MishraCES, 
   vrs = TRUE, rho2 = (0:4)/3-0.6, control = list( maxit = 5000 ) )
print.default( cesBfgsVrsGrid2 )
print( cesBfgsVrsGrid2 )
summary( cesBfgsVrsGrid2 )

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

## Levenberg-Marquardt, Grid Search for rho_1, rho_2, and rho, TC, VRS
cesLmTcVrsGrid123 <- cesEst( "yObs", xNames, data = MishraCES, vrs = TRUE,
   rho1 = (0:2)/2+0.3, rho2 = (0:2)/1.5-0.9, rho = (0:2)/2+0.1,
   control = nls.lm.control( maxiter = 200 ) )
print.default( cesLmTcVrsGrid123 )
print( cesLmTcVrsGrid123 )
summary( cesLmTcVrsGrid123 )

