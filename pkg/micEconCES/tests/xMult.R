# load the micEconCES package
library( micEconCES )

# seed for random number generation
set.seed( 123 )

# number of observations
nObs <- 200

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

# coefficients
cesCoef <- c( 1, 1:nExog / sum( 1:nExog ), 0.5, 1.1 )
names( cesCoef ) <- c( "gamma", paste( "delta", 1:nExog, sep = "_" ),
   "rho", "nu" )

# calculate deterministic endogenous variable
cesData$y <- cesCalc( xNames = xxNames, data = cesData, coef = cesCoef )
print( cesData$y )
# check if removing the names of the coefficients makes a difference
all.equal( cesData$y,
   cesCalc( xNames = xxNames, data = cesData, coef = unname( cesCoef ) ) )
# check if permuting the coefficients makes a difference
all.equal( cesData$y,
   cesCalc( xNames = xxNames, data = cesData, coef = sample( cesCoef, 7 ) ) )

# adding noise to the endogenous variable
cesData$y <- cesData$y + rnorm( nObs )
