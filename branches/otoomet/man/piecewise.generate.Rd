\name{piecewise.generate}
\alias{piecewise.generate}
\title{Generate survival data}
\description{
  Generates two-state single-spell survival data with
  piecewise-constant-baseline hazard.  Intended to be modified to
  satisfy personal needs.
}
\usage{
piecewise.generate(N, het=c("discrete", "gamma", "none"), NX=3,
cens=TRUE, cens.time=c(5, 20))
}
\arguments{
  \item{N}{Number of observations to generate}
  \item{het}{type of unobserved heterogeneity to be generated.
    \code{none} means no heterogeneity}
  \item{NX}{number of observable regressors}
  \item{cens}{logical, whether to censore observations}
  \item{cens.time}{numerical vector of length 2, if \code{cens=TRUE},
    observations are censored on random times, uniformly distributed
    between those values.  If you want a fixed censoring time, use
    \code{cens.time=c{10,10}.}}
}
\details{"Observable" covariates are generated as \code{runif(N)}.
  "Unobservable" covariates are normalised as follows:
  \itemize{
    \item{discrete:} probability 0.3 value \eqn{e^{-1.5}}{exp(-1.5)},
    probability 0.7 value 1.
    \item{gamma:} gamma distribution with mean 1 and variance 0.25.
  }
  Time intervals are generated as 0-2, 2-5 and 5-\ldots
  }
\value{A list with following components:
  \item{X}{\eqn{N\times NX}{NxNX} matrix of observable covariates}
  \item{v}{vector of heterogeneity distributet according to \code{het}.
    This is supposed to be "unobservable"}
  \item{t}{vector of (possibly censored) survival times}
  \item{cens}{logical vector, whether particular observation is censored}
}
\author{Ott Toomet \email{siim@obs.ee}}
\seealso{\code{\link{piecewise.discrete}}}

\examples{
## generate some data with discrete heterogeneity:
h <- piecewise.generate(400)
## And estimate it.  
summary(piecewise.discrete(h$t, h$cens, h$X, c(2,5), iterlim=50))

## Now generate some more, without unobservable heterogeneity
h <- piecewise.generate(1000, het="none", cens.time=c(7,7))
## And estimate it
summary(piecewise(h$t, h$cens, h$X, c(2,5), iterlim=50))
}
\keyword{models}
\keyword{survival}
