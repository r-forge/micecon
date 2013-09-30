library( "micEconNP" )

set.seed( 123 )
n <- 250
x1 <- runif(n)
x2 <- rnorm(n)
z1 <- runif(n)
z2 <- runif(n, min=-2, max=2)
y <- sin( 5 * z1 ) + x1 * exp( z2 ) + x2 * ( z1^2 + z2 ) + rnorm( n, sd = 0.2 )
model <- npscoef( y ~ x1 + x2 | z1 + z2, beta = TRUE )

print( model )
summary( model )
round( coef( model ), 3 )

grad <- npscoefGrad( model )
round( grad, 3 )

gradTrue <- array( 0, c( n, 4, 2 ) )
gradTrue[ , 1, 1 ] <- 5 * cos( 5 * z1 )
gradTrue[ , 3, 1 ] <- 2 * z1
gradTrue[ , 2, 2 ] <- exp( z2 )
gradTrue[ , 3, 2 ] <- 1
zData <- cbind( z1, z2 )
# for( j in 1:2 ) {
#    for( i in 1:3 ) {
#       plot( zData[,j], grad[ , i, j ],
#          main = paste( "d b", i-1, " / d z", j, sep = "" ) )
#       points(  zData[,j], gradTrue[ , i, j ], col = "red" )
#       cat( "Press <enter>\n")
#       readLines( n = 1 )
#    }
# }
all.equal( grad[ , 1, ] + grad[ , 2, ] * x1 + grad[ , 3, ] * x2, grad[ , 4, ] )