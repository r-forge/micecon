\name{NObs}
\alias{NObs}
\alias{NObs.default}
\title{Return number of observations for statistical models} 
\description{
  Returns number of observations for statistical models.  The default
  method assumes presence of a component \code{$param$NObs} in {x}.
}
\usage{
NObs(x, ...)
NObs.defalt(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a statistical model, such as created by \link{\code{lm}}}
  \item{\dots}{further arguments for methods}
}
\details{
  This is a generic function. Functions with names beginning in
  NObs. will be methods for this function.  
}
\value{
  A numeric, number of observations
}
\author{Ott Toomet, \email{otoomet@econ.au.dk}}
\seealso{\code{\link{NParam}}}
\examples{
## Do a trivial probit example and find number of observations
x <- runif(100)
y <- x + rnorm(100) > 0
p <- probit(y ~ x)
NObs(p)
}
\keyword{methods}
