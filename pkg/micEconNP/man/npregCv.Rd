\name{npregCv}
\alias{npregCv}
\alias{npplregCv}

\title{Cross-Validation Nonparametric Models}

\description{
   Leave-one-out cross-validation
   of nonparametric kernel regression models.
}

\usage{
   npregCv( x )

   npplregCv( x )
}

\arguments{
   \item{x}{an object of class \code{npregression}
      (returned by \code{\link[np]{npreg}})
      or an object of class \code{plregression}
      (returned by \code{\link[np]{npplreg}}).}
}


\value{
   A single numeric value:
   the mean squared prediction error over all observations,
   where the predicted values are obtained
   based on models that are estimated without the respective observation.
   This value has an attribute \code{err}
   that is the vector of the prediction errors of all observations.
}

\seealso{\code{\link[np]{npreg}}, \code{\link[np]{npplreg}}.}

\author{Arne Henningsen}

\examples{
set.seed(42)
n <- 250
x1 <- rnorm(n)
x2 <- rnorm(n)
y <- 1 + x1 + x2^2 + rnorm(n)

model <- npreg( y ~ x1 + x2, regtype="ll" )

cv <- npregCv( model )
all.equal( cv, model$bws$fval )

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
cv
}

\keyword{models}
\keyword{nonparametric}
\keyword{regression}
