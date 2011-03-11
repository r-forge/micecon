# load the micEconCES package
library( micEconCES )

# seed for random number generation
set.seed( 123 )

# number of observations
nObs <- 20

# number of explanatory variables
nExog <- 4

# create data set
cesData <- data.frame( obsNo = 1:nObs )

# names of explanatory variables
xxNames <- paste( "xx", 1:nExog, sep = "." )

# add explanatory variables
for( i in 1:nExog ) {
   cesData[[ xxNames[ i ] ]] <- rchisq( nObs, 10 + i )
}
cesData$time <- c( 0:( nObs - 1 ) )

# delta coefficients
cesDelta <- c( 1:nExog / sum( 1:nExog ) )
names( cesDelta ) <- paste( "delta", 1:nExog, sep = "_" )

# vector with values around 0 for coefficient "rho"
rhos = c( -exp(-(1:20)),0,exp(-(20:1)) )

# matrix for returned endogenous variables
y <- matrix( NA, nrow = length( rhos ), ncol = nObs )
rownames( y ) <- c( -(1:20), 0, (20:1) )

# calculate endogenous variables
for( i in 1:length( rhos ) ) {
   # coefficients
   cesCoef <- c( gamma = 1, lambda = 0.02, cesDelta, rho = rhos[ i ], nu = 1.1 )
   y[ i, ] <- cesCalc( xNames = xxNames, tName = "time", data = cesData, 
      coef = cesCoef )
}

# print matrix of endogenous variables
print( y )

# endogenous variables in case of a Cobb-Douglas function (rho = 0)
cdCoef <- c( a_0 = unname( log( cesCoef[ "gamma" ] ) ) )
for( i in 1:nExog ) {
   cdCoef <- c( cdCoef, cesCoef[ paste( "delta", i, sep = "_" ) ] * cesCoef[ "nu" ] )
}
names( cdCoef ) <- paste( "a", 0:nExog, sep = "_" )
yCd <- cobbDouglasCalc( xNames = xxNames, data = cesData, 
   coef = cdCoef ) * exp( cesCoef[ "lambda" ] * cesData$time )

# print endogenous variables for different rhos (adjusted with the y at rho=0)
for( i in 1:nObs ) {
   print( format( round( y[ , i, drop = FALSE ] - yCd[i], 11 ), 
      scientific = FALSE ) )
}
