\name{summary.tobit5}
\alias{summary.tobit5}
\alias{print.summary.tobit5}
\title{Summary and print methods for tobit-models}
\description{
  Summarizing tobit-model fits
}
\usage{
summary.tobit5(object, ...)
print.summary.tobit5(x, ...)
}
\arguments{
  \item{object}{object of class "tobit5", usually result of a call to
    \code{\link{tobit5}}.}
  \item{x}{object of class "summary.tobit5", usually a result of a call
    to \code{summary.lm}.}
  \item{\dots}{currently ignored}
}
\details{
  \code{summary.tobit5} calculates standard errors, t-values and
  corresponding probabilities, based on variance-covariance matrix.
  \code{print.summary.tobit5} prints all those in a nice form.
  }
\value{
  \code{summary.tobit5} return an object of class "summary.tobit5".  It
  is the same as for \code{\link{summary.maximisation}}}
\author{Ott Toomet \email{siim (a) obs ee}}
\seealso{\code{\link{tobit5}}}
\examples{
## uniform observable covariates (note that exclusion restriction is
## fulfilled)! 
z <- runif(200)
x <- runif(200)
## latent selection variable
y1s <- 2*z + rnorm(200)
## observable choice
y1 <- y1s > 0
## latent regressions
y2 <- x + 0.5*rnorm(200) # observable if y1 = 0
y3 <- -x + 0.5*rnorm(200) # observable if y1 = 1
a <- tobit5(y1~z, y2~x, y3~x)
summary(a)
}
\keyword{methods}
\keyword{print}
