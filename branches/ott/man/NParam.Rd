\name{NParam}
\alias{NParam}
\alias{NParam.default}
\alias{NParam.lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Number of model parameters}
\description{
  This function returns the number of model parameters.  The default
  method returns the component \code{x$param$NParam}.
}
\usage{
NParam(x, ...)
NParam.default(x, \dots)
NParam.lm(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\details{
  Not all parameters need to be free.  The number may depend on the
  normalisation. 
}
\arguments{
  \item{x}{a statistical model}
  \item{\dots}{other arguments for methods}
}
\value{
  Number of parameters in the model
}
\author{Ott Toomet, \email{otoomet@econ.au.dk}}
\seealso{\code{\link{NObs}} for number of observations}
\examples{
# Construct a simple OLS regression:
x1 <- runif(100)
x2 <- runif(100)
y <- 3 + 4*x1 + 5*x2 + rnorm(100)
m <- lm(y~x1+x2)  # estimate it
summary(m)
NParam(m) # you get 3
}
\keyword{methods}
