library( "miscTools" )

set.seed( 123 )
m1 <- matrix( rnorm( 9 ), ncol = 3 )
print( m1 )
semidefiniteness( m1 )
semidefiniteness( m1, positive = FALSE )

m2 <- crossprod( m1 )
print( m2 )
semidefiniteness( m2 )
semidefiniteness( m2, positive = FALSE )
semidefiniteness( -m2 )
semidefiniteness( -m2, positive = FALSE )

m3 <- cbind( m2, - rowSums( m2 ) )
m3 <- rbind( m3, - colSums( m3 ) )
print( m3 )
semidefiniteness( m3 )
semidefiniteness( m3, positive = FALSE )

m4 <- m3 * 1e6
print( m4 )
# rcond(m4)
# det(m4)
semidefiniteness( m4 )
semidefiniteness( m4, positive = FALSE )

m5 <- diag( -1, 4, 4 )
semidefiniteness( m5 )
semidefiniteness( m5, positive = FALSE )

m6 <- matrix( -1, 4, 4 )
semidefiniteness( m6 )
semidefiniteness( m6, positive = FALSE )
