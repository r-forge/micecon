\name{activePar}
\alias{activePar}
\alias{activePar.default}
\title{free parameters under maximisation}
\description{
  Return a logical vector, indicating which parameters were free under
  maximisation, as opposed to parameters, fixed by user.
}
\usage{
activePar(x, \dots)
activePar.default(x, \dots)
}
\arguments{
  \item{x}{object, created by a maximisation routine, or derived from a
    maximisation object.  Currently only \code{\link{maxNR}} and it's
    derivations support \code{activePar}}
  \item{\dots}{further arguments for methods}
}
\details{
  Several optimisation routines allow the user to fix some parameter
  values (or do it automatically in some cases).  One has to know know
  which parameters are free for any gradient or Hessian based inference,
  e.g. for inverting Hessian in order to find standard errors.
}
\value{
  A logical vector, indicating whether the parameters were free to
  change during optimisation algorithm.
}
\author{Ott Toomet, \email{otoomet@econ.au.dk}}
\seealso{\code{\link{maxNR}}, \code{\link{NObs}}}
\examples{
# a simple two-dimensional exponential hat
f <- function(a) exp(-a[1]^2 - a[2]^2)
#
# maximize wrt. both parameters 
free <- maxNR(f, theta=1:2) 
summary(free)  # results should be close to (0,0)
# allow only the second parameter to vary
cons <- maxNR(f, theta=1:2, activePar=c(FALSE,TRUE))
summary(cons) # result should be around (1,0)
}
\keyword{methods}
\keyword{optimize}
