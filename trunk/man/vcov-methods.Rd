\name{vcov-methods}
\alias{vcov.probit}
\alias{vcov.heckit}
\alias{vcov.maxLik}

\title{Methods for Function 'vcov' in Package 'micEcon'}
\description{
   Extract variance-covariance matrices.
}
\usage{
\method{vcov}{heckit}(object, part = "outcome", ...)
   \method{vcov}{probit}( object, ... )
   \method{vcov}{maxLik}( object, ... )
}
\arguments{
   \item{object}{an object of class \code{probit}, \code{heckit},
      or \code{maxLik}.}
  \item{part}{character string: which parts to extract: "full" for all
    the estimated parameters (probit selection estimates, outcome
    estimates, including inverse Mills ratio), error variance and
    correlation; or "outcome" for the outcome model results only
    (including the inverse Mills ratio)}
   \item{\dots}{further arguments (currently ignored).}
}
\details{
  variance-covariance matrix of the heckit model is only partially
  implemented.  The rest of
  the matrix is filled with NA-s.
}
\value{
   the estimated variance covariance matrix of the coefficients.
}
\author{
   Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de},
   Ott Toomet \email{otoomet@ut.ee}
}

\seealso{\code{\link[stat]{vcov}}, \code{\link{maxLik}}.}
\examples{
## Estimate a simple female wage model taking into account the labour
## force participation
   data(Mroz87)
   a <- heckit(lfp ~ huswage + kids5 + mtr + fatheduc + educ + city,
               log(wage) ~ educ + city, data=Mroz87)
## print the variance-covariance matrix of the outcome model:
print(vcov(a))
## now print the full matrix (partly unimplemented)
print(vcov(a, part="full"))
}
\keyword{methods}
