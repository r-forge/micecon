\name{compare.derivatives}
\alias{compare.derivatives}
\title{function to compare analytic and numeric derivatives}
\description{
  This function compares analytic and numerical derivative and prints
  diagnostics.  Intended for testing pre-programming derivative routines
  for maximisation algorithm.
}
\usage{
compare.derivatives(f, grad, hess=NULL, t0, eps=1e-6, ...) 
}
\arguments{
  \item{f}{Function to be differentiated.  The parameter of interest
    must be the first argument.  May return a vector.}
  \item{grad}{Analytic gradient.  Must use the same parameters as f.  If
    f is a vector-valued function, grad must return a matrix where the
    number of rows equals the number of components of f, and the number
    of columns must equal to the number of components in t0}
  \item{hess}{Analytic hessian.  If present, hessian matrices are
    compared too.  Only appropriate for scalar-valued function.}
  \item{t0}{Parameter vector.  The derivative must be taken w.r.t. this
    vector.}
  \item{eps}{Interval to be used for numerical derivatives.  Central
    derivative is used.}
  \item{\dots}{further arguments to f, grad and hess.}
}
\details{
  For every component of f, the parameter value, analytic and numeric
  derivative and their relative difference (analytic - numeric)/analytic
  are printed.  If analytic derivatives are correct and the function
  looks nice, expect the relative differences to be less than 1e-7.
}
\author{Ott Toomet \email{siim@obs.ee}}
\seealso{\code{\link{numeric.gradient}}}

\examples{
## A simple example with sin(x)' = cos(x)
f <- sin
compare.derivatives(f, cos, t0=1) 
##
## Example of log-likelihood of normal density.  Two-parameter
## function. 
x <- rnorm(100, 1, 2) # generate rnorm x
l <- function(b) sum(log(dnorm((x-b[1])/b[2])/b[2]))
              # b[1] - mu, b[2] - sigma
gradl <- function(b) {
c(sum(x - b[1])/b[2]^2,
sum((x - b[1])^2/b[2]^3 - 1/b[2]))
}
compare.derivatives(l, gradl, t0=c(1,2))
}
\keyword{math}
\keyword{utilities}
