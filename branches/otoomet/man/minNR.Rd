\name{minNR}
\alias{minNR}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Newton-Raphson Minimization}
\description{Finds a minimum of the given function using Newton-Raphson
  algorithm. 
}
\usage{
minNR(f, theta0, names=NULL, print.level=0, tol=1e-06, gradtol=1e-06,
steptol=1e-06, lambdatol=1e-06, iterlim=15, \dots)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{f}{function to be minimized.  First argument must be the
    parameter vector.  Single number with attributes
    \code{gradient} and \code{hessian}}
  \item{theta0}{vector of initial parameters}
  \item{names}{vector of parameter names.  Used in output routines}
  \item{print.level}{0 - nothing is printed, 1 - initial and final
    values, 2 - parameter values at every iteration, 3 - parameters and
    hessian at every iteration}
  \item{tol}{stop if successive function values differ less than \code{tol}}
  \item{gradtol}{stop if gradient smaller than \code{gradtol}}
  \item{steptol}{when whole step into calculated direction do not result
    in lower function value, half of the previous step is tried.  This
    process is continued until step < steptol.  Then exit code is set to
    3}
  \item{lambdatol}{it is required that the minimal eigenvalue of hessian
    is larger than \code{lambdatol}.  Otherwise, a small diagonal matrix
    is addet to the hessian, so that resulting sum satisfies the
    requirement}
  \item{iterlim}{maximum number of iterations}
  \item{\dots}{additional arguments for \code{f}}
}
\value{List of class "maximization".  There is a \code{\link{summary}}
  method available.  List contains following components:
  \item{maximum}{minimum of the function}
  \item{estimate}{vector of estimated parameter values}
  \item{names}{vector of parameter names}
  \item{gradient}{gradient at minimum (normally close to 0)}
  \item{hessian}{Hessian matrix at minimum}
  \item{code}{return code}
  \item{message}{message, describing the return code}
  \item{last.step}{only if \code{code == 3} (not able to locate a better
    point near current estimation), otherwise \code{NULL}:
    \itemize{
      \item{theta0:}{last parameter estimate}
      \item{f0:}{last function value}
      \item{theta1:}{new parameter estimate}
    }
  }
  \item{iterations}{actual number of iterations}
  \item{type}{type of optimization (= "Newton-Raphson minimization")}
}
\references{Amemiya T. (1985), \emph{"Advanced Econometrics"}, Harvard.}
\author{Ott Toomet \email{siim@obs.ee}}
\note{This function is essentially the same as \code{nlm}}
\section{WARNING}{You should not use numerical gradient and hessian in
  the algorith.  They are much more unstable.}
\seealso{\code{\link{optim}}, \code{\link{nlm}}}

\examples{
## a 2D quadratic function:
f <- function(b) { x <- b[1]; y <- b[2];
    val <- (x - 2)^2 + (y - 3)^2
    attr(val, "gradient") <- c(2*x - 4, 2*y - 6)
    attr(val, "hessian") <- matrix(c(2, 0, 0, 2), 2, 2)
    val
}
## Note that NR finds the minimum of a quadratic function with a single
## iteration.  Use c(0,0) as initial value:
summary(minNR(f, c(0,0)))
## Now use c(1000000, -777777) as initial value:
summary(minNR(f, c(1000000,-777777)))
}
\keyword{optimize}% at least one, from doc/KEYWORDS

