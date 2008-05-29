\name{coef.quadFuncEst}
\alias{coef.quadFuncEst}

\title{Coefficients of a Quadratic Function}

\description{
   This method returns the coefficients from the results
   of an estimated quadratic function.
}

\usage{
   \method{coef}{quadFuncEst}( object, ... )
}

\arguments{
  \item{object}{an object of class \code{quadFuncEst}.}
  \item{...}{currently ignored.}
}

\value{
The \code{coef} method returns a vectors containing all (linearly independent)
coefficients of a quadratic function.}
}

\author{Arne Henningsen, \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{quadFuncEst}}}

\keyword{methods}
