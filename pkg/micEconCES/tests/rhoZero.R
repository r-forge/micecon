# load the micEconCES package
library( micEconCES )

# seed for random number generation
set.seed( 123 )

# number of observations
nObs <- 20

# create data set with explanatory variables
cesData <- data.frame( xx1 = rchisq( nObs, 10 ), xx2 = rchisq( nObs, 10 ) )

# names of explanatory variables
xxNames <- c( "xx1", "xx2" )

# vector with values around 0 for coefficient "rho"
rhos = c( -exp(-(1:20)),0,exp(-(20:1)) )

# matrix for returned endogenous variables
y <- matrix( NA, nrow = length( rhos ), ncol = nObs )
rownames( y ) <- c( -(1:20), 0, (20:1) )

# calculate endogenous variables
for( i in 1:length( rhos ) ) {
   # coefficients
   cesCoef <- c( gamma = 1, delta = 0.6, rho = rhos[ i ], nu = 1.1 )
   y[ i, ] <- cesCalc( xNames = xxNames, data = cesData, coef = cesCoef )
}

# print matrix of endogenous variables
print( y )

# endogenous variables in case of a Cobb-Douglas function (rho = 0)
cdCoef <- c( a_0 = unname( log( cesCoef[ "gamma" ] ) ), 
   a_1 = unname( cesCoef[ "delta" ] * cesCoef[ "nu" ] ),
   a_2 = unname( ( 1 - cesCoef[ "delta" ] ) * cesCoef[ "nu" ] ) )
yCd <- cobbDouglasCalc( xNames = xxNames, data = cesData, 
   coef = cdCoef )

# print endogenous variables for different rhos (adjusted with the y at rho=0)
for( i in 1:nObs ) {
   print( format( round( y[ , i, drop = FALSE ] - yCd[i], 11 ), 
      scientific = FALSE ) )
}
