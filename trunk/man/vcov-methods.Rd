\name{vcov-methods}
\alias{vcov.probit}
\alias{vcov.heckit}
\alias{vcov.maxLik}

\title{Methods for Function 'vcov' in Package 'micEcon'}
\description{
   Extract variance-covariance matrices.
}
\usage{
   \method{vcov}{probit}( object, ... )
   \method{vcov}{heckit}( object, ... )
   \method{vcov}{maxLik}( object, ... )
}
\arguments{
   \item{object}{an object of class \code{probit}, \code{heckit},
      or \code{maxLik}.}
   \item{\dots}{further arguments (currently ignored).}
}
\value{
   the estimated variance covariance matrix of the coefficients.
}
\author{
   Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de},
   Ott Toomet \email{otoomet@ut.ee}
}

\seealso{\code{\link[stat]{vcov}}, \code{\link{maxLik}}.}

\keyword{methods}
