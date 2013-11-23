library( "micEconNP" )
options( np.messages = FALSE )
set.seed(42)
n <- 250
x1 <- rnorm(n)
x2 <- rnorm(n)
y <- 1 + x1 + x2^2 + rnorm(n)

model <- npreg( y ~ x1 + x2, regtype="ll" )

cv <- npregCv( model )
round( c( cv ), 3 )
round( attr( cv, "err" ), 3 )
all.equal( cv, model$bws$fval, check.attributes = FALSE )
