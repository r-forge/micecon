\name{npregCv}
\alias{npregCv}
\alias{npplregCv}
\alias{crsCv}

\title{Cross-Validation Nonparametric Models}

\description{
   Leave-one-out cross-validation
   of nonparametric kernel regression models,
   partial linear models, and categorical regression splines.
}

\usage{
   npregCv( x )

   npplregCv( x )
   
   crsCv( x )
}

\arguments{
   \item{x}{an object of class \code{npregression}
      (returned by \code{\link[np]{npreg}})
      or an object of class \code{plregression}
      (returned by \code{\link[np]{npplreg}})
      or an object of class \code{crs}
      (returned by \code{\link[crs]{crs}}).}
}


\value{
   A single numeric value:
   the mean squared prediction error over all observations,
   where the predicted values are obtained
   based on models that are estimated without the respective observation.
   This value has an attribute \code{err}
   that is the vector of the prediction errors of all observations.
}

\seealso{\code{\link[np]{npreg}}, \code{\link[np]{npplreg}},
   \code{\link[crs]{crs}}.}

\author{Arne Henningsen}

\examples{
set.seed(42)
n <- 250
x1 <- rnorm(n)
x2 <- rnorm(n)
y <- 1 + x1 + x2^2 + rnorm(n)

model <- npreg( y ~ x1 + x2, regtype="ll" )

cv <- npregCv( model )
c( cv )
all.equal( cv, model$bws$fval, check.attributes = FALSE )


# partially linear model
set.seed(42)
n <- 250
x1 <- rnorm(n)
x2 <- rbinom(n, 1, .5)
z1 <- rbinom(n, 1, .5)
z2 <- rnorm(n)
y <- 1 + x1 + x2 + z1 + sin(z2) + rnorm(n)

model <- npplreg( y ~ x1 + factor(x2) | factor(z1) + z2, regtype="ll" )

cv <- npplregCv( model )
c( cv )


# categorical regression splines
set.seed(42)
n <- 250
num.eval <- 50
x1 <- runif(n)
x2 <- runif(n)
z <- round( runif( n, min = 0, max = 3 ) )
dgp <- cos( 2 * pi * x1 ) + sin( 2 * pi * x2 ) + z/5
z <- factor(z)
y <- dgp + rnorm( n, sd = 0.5 )

model <- crs( y ~ x1 + x2 + z, deriv = 1, nmulti = 1 )

cv <- crsCv( model )
c( cv )
all.equal( cv, model$cv.min, check.attributes = FALSE )
}

\keyword{models}
\keyword{nonparametric}
\keyword{regression}
