\name{selection}
\alias{selection}
\title{Heckman-style selection models}
\description{
  This functions calculates Heckman-style selection models either with
  one or two outcomes.
}
\usage{
selection(selection, outcome, data = sys.frame(sys.parent()),
method = "ml", init = NULL, 
ys = FALSE, xs = FALSE,
yo = FALSE, xo = FALSE,
mfs = FALSE, mfo = FALSE,
print.level = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{selection}{formula, the selection equation}
  \item{outcome}{the outcome equation(s).  Either a single equation (for
  tobit 2 models), or a list of two equations (tobit 5 models).}
  \item{data}{an optional data frame, list or environment (or object
          coercible by \code{as.data.frame} to a data frame) containing the
          variables in the model.  If not found in \code{data}, the
          variables are taken from \code{environment(formula)}, typically
          the environment from which \code{selection} is called.}
  \item{method}{how to estimate the model.  Either "ml" for Maximum
    Likelihood, or "2step" for 2-step estimation}
  \item{init}{initial values for the ML estimation.  Should have names}
  \item{ys, yo, xs, xo, mfs, mfo}{logicals.  If true, the response (\code{y}),
    model matrix (\code{x}) or the model frame (\code{mf})
    of the selection (\code{s}) or outcome
    (\code{o}) equation(s) are returned}
  \item{print.level}{integer.  Various debugging information, higher
    value gives more information.}
  \item{\dots}{additional parameters for the corresponding fitting
    functions \code{\link{tobit2fit}} and \code{\link{tobit5fit}}}
}
\details{
  These selection models assume a known (multivariate normal)
  distribution
  of error terms.  Because of this, the instruments (exclusion
  restrictions) are not necessary.  However, if no instruments are
  supplied, the results are based solely on the assumption on
  multivariate normality.  This may or may not be an appropriate
  assumption for a particular problem.
}
\value{
  object of class 'sampleSelection', which inherits from class
  'maxLik'.  A list, which has all the components of 'maxLik', and in
  addition
  \item{twoStep}{(only if initial values not given)
    results of the 2-step estimation, used for initial values}
  \item{init}{initial values for ML estimation}
  \item{param}{various parameters, including ....}
  \item{termsS, termsO}{terms for the selection and outcome equation}
  \item{ys, xs, yo, xo, mfs, mfo}{response, matrix and frame of the
    selection- and outcome equations (as a list of two for the latter).
    NULL, if not requested.}
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
## tobit-5 example
N <- 500
   library(mvtnorm)
   vc <- diag(3)
   vc[lower.tri(vc)] <- c(0.9, 0.5, 0.6)
   vc[upper.tri(vc)] <- vc[lower.tri(vc)]
   eps <- rmvnorm(N, rep(0, 3), vc)
   xs <- runif(N)
   ys <- xs + eps[,1] > 0
   xo1 <- runif(N)
   yo1 <- xo1 + eps[,2]
   xo2 <- runif(N)
   yo2 <- xo2 + eps[,3]
   a <- selection(ys~xs, list(yo1 ~ xo1, yo2 ~ xo2))
summary(a)
##
## tobit2 example
   vc <- diag(2)
   vc[2,1] <- vc[1,2] <- -0.7
   eps <- rmvnorm(N, rep(0, 2), vc)
   xs <- runif(N)
   ys <- xs + eps[,1] > 0
   xo <- runif(N)
   yo <- (xo + eps[,2])*(ys > 0)
   a <- selection(ys~xs, yo ~xo)
summary(a)
}

}
\keyword{models}
\keyword{regression}
