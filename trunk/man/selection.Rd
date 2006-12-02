\name{selection}
\alias{selection}
\title{Heckman-style selection models}
\description{
  This functions calculates Heckman-style selection models either with
  one or two outcomes.
}
\usage{
selection(selection, outcome, data = sys.frame(sys.parent()),
method = "ml", init = NULL, print.level = 0,
ySelection = FALSE, xSelection = FALSE,
yOutcome = FALSE, xOutcome = FALSE,
model = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{selection}{formula, the selection equation}
  \item{outcome}{ ~~Describe \code{outcome} here~~ }
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{method}{ ~~Describe \code{method} here~~ }
  \item{init}{ ~~Describe \code{init} here~~ }
  \item{print.level}{ ~~Describe \code{print.level} here~~ }
  \item{ySelection}{ ~~Describe \code{ySelection} here~~ }
  \item{xSelection}{ ~~Describe \code{xSelection} here~~ }
  \item{yOutcome}{ ~~Describe \code{yOutcome} here~~ }
  \item{xOutcome}{ ~~Describe \code{xOutcome} here~~ }
  \item{model}{ ~~Describe \code{model} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
  These selection models assume known (multivariate normal) distribution
  of error terms.  Because of this. instruments (exclusion
  restrictions) are not necessary.  However, it no instruments are
  supplied, the results are based solely on the assumption on
  multivariate normality.  This may or may not be an appropriate
  assumption for a particular problem.
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
