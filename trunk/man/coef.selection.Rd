\name{coef.selection}
\alias{coef.selection}
\title{Extract Coefficients from Selection Models}
\description{
  This function extracts coefficients
  from sample selection models
}
\usage{
\method{coef}{selection}(object, part = "full", ...)
}
\arguments{
  \item{object}{object of class "selection".}
  \item{part}{character string indication which parts to extract:
   "full" for all estimated parameters (selection estimates,
   outcome estimates, error variance and correlation)
   or "outcome" for the outcome estimates only
   (including the coefficient of the inverse Mill's ratio
   in case of a two-step estimation).}
  \item{\dots}{currently not used.}
}
\value{
  numeric vector of coefficients.
}
\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de},
   Ott Toomet (\email{otoomet@ut.ee})}

\seealso{\code{\link{selection}}, \code{\link{coefficients}},
   and \code{\link{vcov}}}

\examples{
## Estimate a simple female wage model taking into account the labour
## force participation
   data(Mroz87)
   a <- heckit(lfp ~ huswage + kids5 + mtr + fatheduc + educ + city,
               log(wage) ~ educ + city, data=Mroz87)
## extract all coefficients of the model:
coef( a )

## now extract the coefficients of the outcome model only:
coef( a, part="outcome")
}
\keyword{methods}
