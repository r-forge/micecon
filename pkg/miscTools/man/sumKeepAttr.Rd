\name{sumKeepAttr}
\alias{sumKeepAttr}

\title{Sum of an Array While Keeping its Attributes}

\description{
   This function returns the sum of an numeric array (e.g. vector or matrix)
   while keeping its attributes.
}

\usage{
sumKeepAttr( x, na.rm = FALSE )
}

\arguments{
   \item{x}{an numeric array (e.g. vector or matrix).}
   \item{na.rm}{logical. Passed to \code{\link[base]{sum}}.
      Should missing values be removed?}
}

\value{
the sum (see \code{\link[base]{sum}}).
}

\author{Arne Henningsen}

\seealso{\code{\link[base]{sum}}}

\examples{
a <- 1:10
attr( a, "min" ) <- 1
attr( a, "max" ) <- 10
sum(a)
sumKeepAttr(a)
}
\keyword{methods}
