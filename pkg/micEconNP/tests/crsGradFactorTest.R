library( "micEconNP" )

## taken from the documentation of crs(), slightly modified

set.seed(42)
## Example - simulated data
n <- 1000
num.eval <- 50
x1 <- runif(n)
x2 <- runif(n)
z <- round( runif( n, min = 0, max = 3 ) )
dgp <- cos( 2 * pi * x1 ) + sin( 2 * pi * x2 ) + z
z <- factor(z)
y <- dgp + rnorm( n, sd = 0.5 )

## Estimate a model with specified degree and bandwidth
model <- crs( y ~ x1 + x2 + z, degree = c(5,5), lambda = c(0.1),
   cv = "none", kernel = TRUE, deriv = 1 )
summary( model )

zGrad <- crsGradFactor( model, "z" )

round( zGrad, 2 )

all.equal( rowSums( zGrad, na.rm = TRUE )[ z != 0 ], 
   model$deriv.mat[ z != 0 , 3 ] )

