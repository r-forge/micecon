\name{numericGradient}
\alias{numericGradient}
\alias{numericHessian}

\title{Functions to Calculate Numeric Derivatives}
\description{
  Calculate (central) numeric gradient and hessian of a function.
  \code{numericGradient} accepts vector-valued functions.
}
\usage{
numericGradient(f, t0, eps=1e-06, ...)
numericHessian(f, grad=NULL, t0, eps=1e-06, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{f}{function to be differentiated.  The first argument must be
    the parameter vector with respect to which it is differentiated.}
  \item{grad}{function, gradient of f}
  \item{t0}{vector, the value of parameters}
  \item{eps}{numeric, the step for numeric differentiation}
  \item{\dots}{furter arguments for \code{f}}
}
\details{
  If the argument \code{grad} is present, hessian takes the numeric
  gradient of the gradient, otherwise the numeric second derivative of
  the \code{f}.
}
\value{
  Matrix.  For \code{numericGradient}, the number of rows is equal to the
  length of the function value vector, and the number of columns is
  equal to the length of the parameter vector.

  For the \code{numericHessian}, both numer of rows and columns is
  equal to the length of the parameter vector.
}
\section{Warning}{
  You should not use numerical differentiation in optimisation
  routines.  Although quite precise in simple cases, they may work really
  poorly in difficult situations.
  }
\author{Ott Toomet \email{siim@obs.ee}}
\seealso{\code{\link{compare.derivatives}}, \code{\link{deriv}}}

\examples{
# A simple example with Gaussian bell
f0 <- function(t0) exp(-t0[1]^2 - t0[2]^2)
numericGradient(f0, c(1,2))
numericHessian(f0, t0=c(1,2))

# An example with the analytic gradient
gradf0 <- function(t0) -2*t0*f0(t0)
numericHessian(f0, gradf0, t0=c(1,2))
# The results should be similar as in the previous case

# The central numeric derivatives have usually quite a high precision
compare.derivatives(f0, gradf0, t0=1:2)
# The differenc is around 1e-10
}
\keyword{math}
\keyword{utilities}
