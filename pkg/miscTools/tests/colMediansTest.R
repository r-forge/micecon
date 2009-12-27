library( "miscTools" )


## matrix
m <- matrix( 1:24, nrow = 6, ncol = 4 )

cm1 <- colMedians( m )
print( cm1 )

rm1 <- rowMedians( m )
print( rm1 )

all.equal( cm1, rowMedians( t( m ) ) )
all.equal( rm1, colMedians( t( m ) ) )


## data.frame
data( "Electricity", package = "Ecdat" )
Electricity <- Electricity[ 1:20, ]

cm2 <- colMedians( Electricity )
print( cm2 )

rm2 <- rowMedians( Electricity )
print( rm2 )

all.equal( cm2, rowMedians( t( Electricity ) ) )
all.equal( rm2, colMedians( t( Electricity ) ) )
