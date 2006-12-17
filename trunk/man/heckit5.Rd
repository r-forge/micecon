\name{heckit5}
\alias{heckit5}
\title{Heckman two-step estimator for tobit5 type of models}
\description{
  This function estimates tobit5 type model by Heckman two-step method.
  It is intended to be called by \code{\link{selection}} for calculating
  the initial values for ML estimation.
}
\usage{
heckit5(selection, outcome1, outcome2, data = sys.frame(sys.parent()),
ySelection = FALSE, xSelection = FALSE, yOutcome = FALSE, xOutcome =
FALSE,
model = FALSE, print.level = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{selection}{formula, the selection equation.  For details, see \code{\link{formula}}.}
  \item{outcome1}{formula, the first outcome equation}
  \item{outcome2}{formula, the second outcome equation}
  \item{data}{an optional data frame containing the variables in the model. If not found in data, the
   variables are taken from environment(formula), typically the
   environment from which lm is called.}
  \item{ySelection}{logical, whether to return the response of the
    selection equation}
  \item{xSelection}{logical, whether to return the model frame for the
    selection equation}
  \item{yOutcome}{logical, whether to return the response of the
    outcome equation}
  \item{xOutcome}{logical, whether to return the model frame for the
    outcome equation}
  \item{model}{ ~~Describe \code{model} here~~ }
  \item{print.level}{integer, larger number gives more debugging information}
}
\details{
  The estimation is proceeded in two steps, analogous to Heckman (1976):
  first, the selection model is estimated by \code{\link{probit}}.
  Thereafter, the outcome equation are estimated, where the expected
  value of the error term is added as an additional explanatory variable.}
\value{
  Object of class "heckit5":
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
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

