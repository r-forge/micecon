\name{summary.cesEst}
\alias{summary.cesEst}
\alias{print.summary.cesEst}
\title{Summarize Estimation of a CES Function}

\description{
   \code{summary} method for objects of class \code{\link{cesEst}}.
}

\usage{
\method{summary}{cesEst}( object, \dots )
\method{print}{summary.cesEst}( x, digits = max(3, getOption("digits") - 3),
   \dots )
}

\arguments{
   \item{object}{an objects returned by \code{\link{cesEst}}.}
   \item{x}{an objects returned by \code{summary.cesEst}.}
   \item{digits}{number of digits.}
   \item{\dots}{further arguments are currently ignored.}
}

\value{
   \code{summary.cesEst} returns a list of class \code{summary.cesEst}
   that contains the elements of the provided \code{object} with with
   following changes or additions:
   \item{coefficients}{a matrix with four columns:
      the estimated coefficients/parameters of the CES
      (including a possible fixed \eqn{\rho}{rho}),
      their standard errors, the t-statistic,
      and corresponding (two-sided) P-values.}
   \item{sigma}{square root of the estimated (asymptotic) variance
      of the random error
      (calculated without correcting for degrees of freedom).}
   \item{r.squared}{\eqn{R^2}-value, i.e. the
      \sQuote{fraction of variance explained by the model}.
      If the model was estimated 
      with argument \code{multErr} equal to \code{TRUE},
      the \eqn{R^2}-value measures the fraction of the explained
      variance of the \emph{logarithmized} dependent variable.}
   \item{vcov}{covariance matrix of the estimated parameters
      (including a possible fixed \eqn{\rho}{rho}).}
}

\seealso{\code{\link{cesEst}} and \code{\link{cesCalc}}.}

\author{Arne Henningsen}

\examples{
   data( germanFarms )
   # output quantity:
   germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
   # quantity of intermediate inputs
   germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput


   ## CES: Land & Labor
   cesLandLabor <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms )

   # print summary results
   summary( cesLandLabor )
}

\keyword{models}
