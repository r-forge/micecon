\name{piecewise}
\alias{piecewise.discrete}
\alias{piecewise}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Duration Models with Piecewise Discrete Baseline Hazard}
\description{
  Two-state single spell duration model with
  proportional hazard and
  piecewise-constant-baseline hazard, without and with two-point
  discrete unobserved heterogeneity.}
}
\usage{
piecewise(t, censored, X=NULL, cm, b0=NULL, \dots)
piecewise.discrete(t, censored, X=NULL, cm, b0=NULL, normalise=FALSE, \dots)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{t}{time, spent in initial state (numeric vector)}
  \item{censored}{which observations are censored (logical vector)}
  \item{X}{matrix of regressors or \code{NULL} if there are no regressors}
  \item{cm}{numeric vector, giving the boundaries of time intervals.
    Note that there will always be one more time interval than the
    number of boundaries, the first starts from 0 and the last ends at
    infinity.}
  \item{b0}{vector for initial parameters.  The length of the parameter
    vector must be \code{3 + ncol(X) + length(cm)}}
  \item{normalise}{whether to normalise the heterogeneity term so that
    it is always less than 1.  One value of heterogeneity term \eqn{v}
    is normalised to be 1,
    the other may be smaller or larger than 1, depending on the data and
    initial parameters (default algorithm gives usually values less than
    1).  This parameter forces that to be smaller than 1, renormalising
    baseline hazard accordingly.
  }
  \item{\dots}{additional argument to \code{minNR} optimization algorithm}
}
\details{
  It is assumed that flow sampling is used to get the data (see
  Lancaster 1990).  The model is estimated with Newton-Raphson
  minimization algorithm.

  The individual hazard is assumed to be in a conventional form:
  \deqn{\theta(t) = v_i e^{\lambda_s} e^{\beta' x_i}}{%
    theta(t) = v[i] exp(lambda[s]) exp(beta' x[i])}
  where \eqn{v_i}{v[i]} is individual unobserved heterogeneity term (not
  present in \code{piecewise()}, \eqn{x_i}{x[i]}
  are covariates and \eqn{\lambda_s}{lambda[s]} describes the
  baseline-hazard for time interval \eqn{s}.  Estimation gives \eqn{\log
    v}{log v} and corresponding logistic probability for one of two
  possible \eqn{v} values.  The other is normalised to 1.
}
\value{List of class "maximization".  There is a
  \code{\link{summary}} method available.  The list has following
  components: 
  \item{maximum}{maximum of the log-likelihood function}
  \item{estimate}{vector of estimated parameter values}
  \item{names}{vector of parameter names}
  \item{gradient}{gradient at maximum (normally close to 0)}
  \item{hessian}{Hessian matrix at maximum}
  \item{code}{return code of maximization algorithm}
  \item{message}{message, describing the return code}
  \item{last.step}{only if \code{code == 3} (not able to locate a better
    point near current estimation), otherwise \code{NULL}:
    \itemize{
      \item{theta0:}{last parameter estimate}
      \item{f0:}{last function value}
      \item{theta1:}{new parameter estimate}
    }
  }
  \item{iterations}{number of iterations}
  \item{type}{test, type of optimization}
}
\references{Lancaster, T. (1990) \emph{"Econometrics of Transition
    Data"},
    Cambridge}
\author{Ott Toomet \email{siim@obs.ee}}
\note{The models is sensitive to censoring.  

  When \code{normalise =
    TRUE}, and the estimate of \eqn{v} was biger than 1, a new parameter
  vector is constructed from the current estimate and one more iteration
  is run.  The new values are used for \code{maximum}, \code{estimate},
  \code{gradient} and \code{hessian}.  Other values (describing the
  maximization process) remain the same.}

\seealso{\code{\link{piecewise.generate}}, \code{\link{minNR}},
  \code{\link{summary.maximization}}}

\examples{
## Generate some data, without unobservable heterogeneity
h <- piecewise.generate(1000, het="none", cens.time=c(7,7))
## And estimate it
summary(piecewise(h$t, h$cens, h$X, c(2,5), iterlim=50))

## Generate some more with discrete heterogeneity:
h <- piecewise.generate(400)
## And estimate it, first as no-heterogeneity case
summary(piecewise(h$t, h$cens, h$X, c(2,5), iterlim=50))
## And now using heterogeneity model.  Note that it takes more iterations:
summary(piecewise.discrete(h$t, h$cens, h$X, c(2,5), iterlim=50))
}
\keyword{models}% at least one, from doc/KEYWORDS
\keyword{survival}
