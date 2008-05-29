\name{coef.translogEst}
\alias{coef.translogEst}

\title{Coefficients of a Translog Function}

\description{
   This method returns the coefficients from the results
   of an estimated translog function.
}

\usage{
   \method{coef}{translogEst}( object, ... )
}

\arguments{
  \item{object}{an object of class \code{translogEst}.}
  \item{...}{currently ignored.}
}

\value{
The \code{coef} method returns a vectors containing all (linearly independent)
coefficients of a translog function.}
}

\author{Arne Henningsen, \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{translogEst}}}

\keyword{methods}
