library( "micEconNP" )
options( np.messages = FALSE )
set.seed(42)
n <- 250
x1 <- rnorm(n)
x2 <- rbinom(n, 1, .5)
z1 <- rbinom(n, 1, .5)
z2 <- rnorm(n)
y <- 1 + x1 + x2 + z1 + sin(z2) + rnorm(n)

model <- npplreg( y ~ x1 + factor(x2) | factor(z1) + z2, regtype="ll" )

cv <- npplregCv( model )
round( c( cv ), 3 )
round( attr( cv, "err" ), 3 )
