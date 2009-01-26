# load micEcon package
library( "micEcon" )

# some tests with the "germanFarms" data set
data( "germanFarms" )

# using all variables
outGF <- outlierHadi( germanFarms[ , -1 ] )
print( outGF )

# only inputs and outputs
outGF2 <- outlierHadi( germanFarms[ ,
   c( "vOutput", "vVarInput", "qLabor", "land" ) ] )
print( outGF2 )
