\name{heckit5}
\alias{heckit5}
\title{Heckman two-step estimator for tobit5 type of models}
\description{
  This function estimates tobit5 type model by the Heckman two-step
  method.
  It is intended to be called by \code{\link{selection}} for calculating
  the initial values for ML estimation.
}
\usage{
heckit5(selection, outcome1, outcome2, data = sys.frame(sys.parent()),
ys = FALSE, yo = FALSE, xs = FALSE, xo = FALSE, mfs = FALSE, mfo = FALSE,
print.level = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{selection}{formula, the selection equation.  For details, see \code{\link{formula}}.}
  \item{outcome1}{formula, the first outcome equation}
  \item{outcome2}{formula, the second outcome equation}
  \item{data}{an optional data frame containing the variables in the model. If not found in data, the
   variables are taken from environment(formula), typically the
   environment from which lm is called.}
  \item{ys, yo, xs, xo, mfs, mfo}{logicals.  If true, the response (\code{y}),
    model matrix (\code{x}) or the model frame (\code{mf})
    of the selection (\code{s}) or outcome
    (\code{o}) equation(s) are returned}
  \item{print.level}{integer, larger number gives more debugging information}
}
\details{
  The estimation is proceeded in two steps, analogous to Heckman (1976):
  first, the selection model is estimated by \code{\link{probit}}.
  Thereafter, the outcome equation are estimated, where the expected
  value of the error term is added as an additional explanatory variable.}
\value{
  Object of class "heckit5": a list with following components:
  \item{probit}{obejct of class "probit", first-step selection estimates}
  \item{twoStep1, twoStep2}{object of class "lm", second-step estimation of the
    outcome equations}
  \item{rho1, rho2}{real, correlation coefficients between the error
    term of the selection equation and the outcome equations}
  \item{sigma1, sigma2}{real, error variance of the outcome equations}
  \item{termsS, termsO}{terms for the selection and outcome equation}
  \item{ys, xs, yo, xo, mfs, mfo}{response, matrix and frame of the
    selection- and outcome equations (as a list of two for the latter).
    NULL, if not requested.}
}
\references{J. Heckman (1976), The common structur of statistical models of
                  truncation, sample selection and limited dependent
                  variables and a simple estimator for such model,
		  Annals of Economic and Social Measurement 5(4), 475-492}
\author{Ott Toomet \email{otoomet@ut.ee}}
\seealso{\code{\link{heckit}}, \code{\link{selection}}}
\examples{
}
\keyword{models}
\keyword{regression}

