\name{summary.maximisation}
\alias{summary.maximisation}
\title{Summary method for maximisation/minimisation}
\description{
  Selects the appropriate information from maximisation algorithm,
  calculates standard errors and
  \eqn{t}-values, and presents in a suitable form.
}
\usage{
summary.maximisation(result, hessian=FALSE, unsucc.step=FALSE)
}
\arguments{
  \item{result}{optimisation result, object of class
    \code{maximisation}.  See \code{\link{minNR}}.}
  \item{hessian}{logical, whether to display Hessian matrix.}
  \item{unsucc.step}{logical, whether to describe last unsuccesful step
    if \code{code} == 3}
}
\value{
  Object of class \code{summary.maximisation}, intended to print with
  corresponding print method.  There are following components:
  \item{type}{type of maximisation.}
  \item{iterations}{number of iterations.}
  \item{code}{exit code (see \code{\link{minNR}}.)}
  \item{message}{a brief message, explaining code.}
  \item{unsucc.step}{description of last unsuccessful step, only if
    requested and \code{code} == 3}
  \item{estimate}{list with following components:
    \itemize{
      \item{value:} function value at maximum
      \item{results:} matrix of results, including estimated values,
      std-deviations, t-values and corresponding probabilities
      \item{hessian:} hessian matrix at maximum, only if requested.
    }
  }
}
\author{Ott Toomet \email{siim@obs.ee}}
\section{WARNING }{Probabilities of t-values depend on several
  assumptions which may not be appropriate in the particular case!}
\seealso{\code{\link{minNR}}}
\examples{
## minimize a 2D quadratic function:
f <- function(b) { x <- b[1]; y <- b[2];
    val <- (x - 2)^2 + (y - 3)^2
    attr(val, "gradient") <- c(2*x - 4, 2*y - 6)
    attr(val, "hessian") <- matrix(c(2, 0, 0, 2), 2, 2)
    val
}
## Note that NR finds the minimum of a quadratic function with a single
## iteration.  Use c(0,0) as initial value.  Note that standard
## deviations, t-values and probabilities are not appropriate here.
summary(minNR(f, c(0,0)))
## Now use c(1000000, -777777) as initial value and ask for hessian
summary(minNR(f, c(1000000,-777777)), hessian=TRUE)
}
\keyword{methods}% at least one, from doc/KEYWORDS
\keyword{print}
