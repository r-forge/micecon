\name{tobit5}
\alias{tobit5}
\title{function to calculate Heckman-type selection models}
\description{
  Heckman-type selection models are regressions with at least two latent
  equation and multivariate normal error term.  See \file{doc/econ.pdf}.
  Tobit-5 model is one particular type of those models.
}
\usage{
tobit5(selection, reg2, reg3, data = sys.frame(sys.parent()),
beta0 = NULL, print.level = 0, ...)
}
\arguments{
  \item{selection}{a symbolic formula for the selection equation}
  \item{reg2}{a symbolic formula for selection 2}
  \item{reg3}{a symbolic formula for selection 3}
  \item{data}{an optional data frame containing the variables in
    the model. By default the variables are
    taken from environment(formula), typically the environment from
    which the function is called}
  \item{beta0}{numeric, initial values of the parameters.  The order is
    as follows: gamma, beta2, sigma2, rho2, beta3, sigma3, rho3 (see
    \file{doc/econ.pdf}).}
  \item{print.level}{information about calculations.  0 -- nothing is
    printed, 1 -- some summary information is printed.
    \code{print.level} is sent further to the maximisation algorithm,
    see \code{\link{maxNR}}.}
  \item{\dots}{further arguments to the maximisation algorithm,
    see \code{\link{maxNR}}.}
}
\details{
  Tobit 5 models have three latent variable regressions with jointly
  normal errors.  If selection variable \eqn{y_1 = 0}{y1 = 0}, selection
  2 is observable, if \eqn{y_1 = 1}{y1 = 1}, selection 3 is observable.
  Values for unobservable selection are undefined, they are commonly
  assumed to be equal to 0.  In the current function the values must be
  there but they are not used.
}
\value{
  a list of class \code{tobit5}:
  \item{results}{results of maximisation of class \code{maximisation}.
    See \code{\link{maxNR}}.}
  \item{probit}{results of the probit-estimation of the selection part
    of the model.  See \code{\link{probit}}.}
}
\references{Amemiya 1985}
\author{Ott Toomet \email{siim@obs.ee}}
\section{WARNING}{These model rely heavily on the parametric
  assumptions.  Do not use them unless the exclusion restriction holds!}
\seealso{\code{\link{probit}} for probit models}
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
# a <- tobit5(y1~z, y2~x, y3~x) # disabled, because tobit5 doesn't 
# summary(a)                    #  work currently.
}
\keyword{models}
\keyword{regression}
