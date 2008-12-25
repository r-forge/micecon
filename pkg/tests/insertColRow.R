## load micEcon package
library( micEcon )

## create a matrix
m <- matrix( 1:9, 3 )

# insert rows
print( insertRow( m, 1, 10:12 ) )
print( insertRow( m, 2, 10:12 ) )
print( insertRow( m, 3, 10:12 ) )
print( insertRow( m, 4, 10:12 ) )

# insert columns
print( insertCol( m, 1, 10:12 ) )
print( insertCol( m, 2, 10:12 ) )
print( insertCol( m, 3, 10:12 ) )
print( insertCol( m, 4, 10:12 ) )

## add row names and column names
rownames( m ) <- c( "R1", "R2", "R3" )
colnames( m ) <- c( "C1", "C2", "C3" )

# insert rows
print( insertRow( m, 1, 10:12 ) )
print( insertRow( m, 2, 10:12 ) )
print( insertRow( m, 3, 10:12 ) )
print( insertRow( m, 4, 10:12 ) )

# insert columns
print( insertCol( m, 1, 10:12 ) )
print( insertCol( m, 2, 10:12 ) )
print( insertCol( m, 3, 10:12 ) )
print( insertCol( m, 4, 10:12 ) )
