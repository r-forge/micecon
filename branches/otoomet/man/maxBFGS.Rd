\name{maxBFGS}
\alias{maxBFGS}

\title{BFGS maximisation}
\description{
  This function is a wrapper for \code{\link{optim}} where the arguments are
  compatible with \code{\link{maxNR}}
}
\usage{
maxBFGS(fn, grad = NULL, theta, print.level = 0, iterlim = 200, tol = 1e-06, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fn}{function to be maximised.  In order to use numeric gradient
    and BHHH method, \code{fn} must return vector of
    observation-specific likelihood values.  Those are summed by maxNR
    if necessary.  If the parameters are out of range, \code{fn} should
    return \code{NA}.  See details for constant parameters.}
  \item{grad}{gradient of the function.  If \code{NULL}, numeric
    gradient is used.  For BHHH method it must return a matrix, where
    rows corresponds to the gradients of the observations.  Note that
    this corresponds to
    \code{t(numericGradient(fn))}, not \code{numericGradient(fn)}.
    It is summed over
    observations in order to get a single gradient vector.}
  \item{theta}{Initial values for the parameters to be optimized over.}
  \item{print.level}{A larger number prints more working information}
  \item{iterlim}{Maximum number of iterations}
  \item{tol}{The absolute convergence tolerance. See \code{\link{optim}}}
  \item{\dots}{Further arguments for \code{fn} and \code{grad}}
}
\value{
  Object of class "maximisation":
  \item{maximum}{The value of fn at maximum}
  \item{estimate}{The best set of parameters found}
  \item{gradient}{The value of gradient at parameter value
    \code{estimate}}
  \item{hessian}{The value of Hessian at optimum}
  \item{code}{Integer success code, 0 is success (see
    \code{\link{optim}})}
  \item{message}{A character string giving any additional information
    returned by the optimizer, or NULL}
  \item{iterations}{A two-element integer vector giving the number of
    calls to fn and gr respectively. This excludes those calls needed to
    compute the Hessian, if requested, and any calls to fn to compute a
    finite-difference approximation to the gradient}
  \item{type}{character string "BFGS maximisation"}
}
\author{Ott Toomet \email{otoomet@ut.ee}}
\seealso{\code{\link{optim}}, \code{\link{nlm}}, \code{\link{maxNR}}, \code{\link{maxBHHH}}}
\examples{
# Maximum Likelihood estimation of Poissonian distribution
n <- rpois(100, 3)
loglik <- function(l) n*log(l) - l - lfactorial(n)
# we use numeric gradient
summary(maxBFGS(loglik, theta=1))
# you would probably prefer mean(n) instead of that ;-)
# Note also that maxLik is better suited for Maximum Likelihood
}
\keyword{optimize}
