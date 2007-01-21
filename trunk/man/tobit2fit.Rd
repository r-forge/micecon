\name{tobit2fit}
\alias{tobit2fit}
\alias{tobit5fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting parametric sample selection models}
\description{
  These models do the actual fitting of the sample selection models.
  The arguments must be given as numeric vectors/matrices, initial value
  of parameters must be specified.  Intended for \code{micEcon} internal use.  
}
\usage{
tobit2fit(YS, XS, YO, XO, init, print.level = 0, ...)
tobit5fit(YS, XS, YO1, XO1, YO2, XO2, init, print.level=0, ...) 
}
\arguments{
  \item{YS}{numeric 0/1 vector, where 0 denotes unobserved outcome
(tobit 2) or outcome 1 observed (tobit 5)}
  \item{XS, XO, XO1, XO2}{numeric matrix, model matrix for selection and
outcome equations}
  \item{YO}{numeric vector, observed outcomes.  Values for unobserved
outcomes are ignored (they may or may not be NA)}
  \item{init}{numeric vector of initial values.  The order is: betaS,
betaO(1), sigma(1), rho(1), betaO2, sigma2, rho2}
  \item{print.level}{numeric, values greater than 0 will produce
increasingly more debugging information}
  \item{\dots}{Additional parameters to \code{\link{maxLik}}}
}
\details{
These functions estimate the models by ML method.
}
\value{
the returned object is identical to \code{\link{maxLik}} returned object
with an extra class \code{tobit2} or \code{tobit5}
}
\author{Ott Toomet \email{otoomet@ut.ee}}
\seealso{\code{\link{selection}}, \code{\link{maxLik}}}
\keyword{internal}
\keyword{models}
\keyword{regression}
