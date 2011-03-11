# load the micEconCES package
library( micEconCES )

# seed for random number generation
set.seed( 123 )

# number of observations
nObs <- 20

# create data set with explanatory variables
cesData <- data.frame( xx1 = rchisq( nObs, 10 ), xx2 = rchisq( nObs, 10 ),
   time = c( 0:( nObs - 1 ) ) )

# names of explanatory variables
xxNames <- c( "xx1", "xx2" )

# vector with values around 0 for coefficient "rho"
rhos <- c( -exp(-(1:20)),0,exp(-(20:1)) )
# rhos <- c( -2^(-(1:40)),0,2^(-(40:1)) )

# matrix for returned endogenous variables
y <- matrix( NA, nrow = length( rhos ), ncol = nObs )
rownames( y ) <- c( -(1:20), 0, (20:1) )

# calculate endogenous variables
for( i in 1:length( rhos ) ) {
   # coefficients
   cesCoef <- c( gamma = 1, lambda = 0.02, delta = 0.6, rho = rhos[ i ], nu = 1.1 )
   y[ i, ] <- cesCalc( xNames = xxNames, tName = "time", data = cesData, 
      coef = cesCoef )
}

# print matrix of endogenous variables
print( y )

# endogenous variables in case of a Cobb-Douglas function (rho = 0)
cdCoef <- c( a_0 = unname( log( cesCoef[ "gamma" ] ) ), 
   a_1 = unname( cesCoef[ "delta" ] * cesCoef[ "nu" ] ),
   a_2 = unname( ( 1 - cesCoef[ "delta" ] ) * cesCoef[ "nu" ] ) )
yCd <- cobbDouglasCalc( xNames = xxNames, data = cesData, 
   coef = cdCoef ) * exp( 0.02 * cesData[[ "time" ]] )

# print endogenous variables for different rhos (adjusted with the y at rho=0)
for( i in 1:nObs ) {
   print( format( round( y[ , i, drop = FALSE ] - yCd[i], 11 ), 
      scientific = FALSE ) )
}


## checking derivatives of the CES with respect to coefficients
# array for returned endogenous variables
deriv <- array( NA, c( length( rhos ), nObs, length( cesCoef ) ) )
dimnames( deriv ) <- list( rhos, 1:nObs, names( cesCoef ) )

# calculate the derivatives
for( i in 1:length( rhos ) ) {
   # coefficients
   cesCoef <- c( gamma = 1, lambda = 0.02, delta = 0.6, rho = rhos[ i ], nu = 1.1 )
   deriv[ i, , ] <- micEconCES:::cesDerivCoef( par = cesCoef, 
      xNames = xxNames, tName = "time", data = cesData, vrs = TRUE, 
      rhoApprox = c( gamma = 5e-6, delta = 5e-6, rho = 1e-3, nu = 5e-6 ) )
}

# print array of derivatives
print( deriv )

# derivatives in case of a Cobb-Douglas function (rho = 0)
derivCd <- matrix( NA, nrow = nObs, ncol = length( cesCoef ) )
dimnames( derivCd ) <- list( 1:nObs, names( cesCoef ) )
# derivCd[ , "gamma" ] <- exp( cesCoef[ "nu" ] *
#    ( cesCoef[ "delta" ] * log( cesData[[ xxNames[ 1 ] ]] ) +
#    ( 1 - cesCoef[ "delta" ] ) * log( cesData[[ xxNames[ 2 ] ]] ) ) )
derivCd[ , "gamma" ] <- 
   cesData[[ xxNames[ 1 ] ]]^( cesCoef[ "nu" ] * cesCoef[ "delta" ] ) *
   cesData[[ xxNames[ 2 ] ]]^( cesCoef[ "nu" ] * ( 1 - cesCoef[ "delta" ] ) ) *
   exp( cesCoef[ "lambda" ] * cesData[[ "time" ]] )
derivCd[ , "lambda" ] <- derivCd[ , "gamma" ] * 
   cesCoef[ "gamma" ] * cesData[[ "time" ]]
derivCd[ , "delta" ] <- cesCoef[ "gamma" ] * cesCoef[ "nu" ] * 
   ( log( cesData[[ xxNames[ 1 ] ]] ) - log( cesData[[ xxNames[ 2 ] ]] ) ) *
   cesData[[ xxNames[ 1 ] ]]^( cesCoef[ "nu" ] * cesCoef[ "delta" ] ) *
   cesData[[ xxNames[ 2 ] ]]^( cesCoef[ "nu" ] * ( 1 - cesCoef[ "delta" ] ) ) *
   exp( cesCoef[ "lambda" ] * cesData[[ "time" ]] )
derivCd[ , "nu" ] <- cesCoef[ "gamma" ] * 
   ( cesCoef[ "delta" ] * log( cesData[[ xxNames[ 1 ] ]] ) +
   ( 1 - cesCoef[ "delta" ] ) * log( cesData[[ xxNames[ 2 ] ]] ) ) *
   cesData[[ xxNames[ 1 ] ]]^( cesCoef[ "nu" ] * cesCoef[ "delta" ] ) *
   cesData[[ xxNames[ 2 ] ]]^( cesCoef[ "nu" ] * ( 1 - cesCoef[ "delta" ] ) ) *
   exp( cesCoef[ "lambda" ] * cesData[[ "time" ]] )
derivCd[ , "rho" ] <- - 0.5 * cesCoef[ "gamma" ] * cesCoef[ "nu" ] *
   cesCoef[ "delta" ] * ( 1 - cesCoef[ "delta" ] ) *
   cesData[[ xxNames[ 1 ] ]]^( cesCoef[ "nu" ] * cesCoef[ "delta" ] ) *
   cesData[[ xxNames[ 2 ] ]]^( cesCoef[ "nu" ] * ( 1 - cesCoef[ "delta" ] ) ) *
   ( log( cesData[[ xxNames[ 1 ] ]] ) - log( cesData[[ xxNames[ 2 ] ]] ) )^2 *
   exp( cesCoef[ "lambda" ] * cesData[[ "time" ]] )



# print derivatives for different rhos (adjusted with the derivatives at rho=0)
for( k in 1:ncol( derivCd ) ) {
   for( i in 1:nObs ) {
      print( format( round( deriv[ , i, k, drop = FALSE ] - 
	 derivCd[ i, k ], 11 ), scientific = FALSE ) )
   }
}
