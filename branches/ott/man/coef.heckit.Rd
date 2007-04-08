\name{coef.heckit}
\alias{coef.heckit}
\alias{coef.heckit5}
\title{results of the heckit model}
\description{
  This functions extract coefficients 
  of the heckit sample selection models
}
\usage{
\method{coef}{heckit}(object, part = "outcome", ...)
\method{coef}{heckit5}(object, part = "outcome", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{object of class "heckit", heckit model}
  \item{part}{character string: which parts to extract: "full" for all
    the estimated parameters (probit selection estimates, outcome
    estimates), error variance and
    correlation; or "outcome" for the outcome model results only}
  \item{\dots}{further arguments for other methods}
}
\value{
  numeric vector of coefficients.
}
\author{Ott Toomet (\email{otoomet@ut.ee})}
\seealso{\code{\link{heckit}}, \code{\link{coefficients}},
  \code{\link{vcov}}, \code{\link{selection}}}
\examples{
## Estimate a simple female wage model taking into account the labour
## force participation
   data(Mroz87)
   a <- heckit(lfp ~ huswage + kids5 + mtr + fatheduc + educ + city,
               log(wage) ~ educ + city, data=Mroz87)
## print the coefficients of the outcome model:
print(coef(a))
## now print all the estimates
print(coef(a, part="full"))
}
\keyword{methods}
