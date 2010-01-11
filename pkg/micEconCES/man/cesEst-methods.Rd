\name{cesEst-methods}
\alias{coef.cesEst}
\alias{coef.summary.cesEst}
\alias{vcov.cesEst}
\title{Methods for Estimated CES Functions}

\description{
Methods for Objects of Class \code{cesEst} and \code{cesEst}.
}

\usage{
\method{coef}{cesEst}( object, \dots )
\method{coef}{summary.cesEst}( object, \dots )
\method{vcov}{cesEst}( object, \dots )
}

\arguments{
   \item{object}{an object of class \code{\link{cesEst}}
      or \code{\link{summary.cesEst}}.}
   \item{\dots}{further arguments are currently ignored.}
}

\value{
   \code{coef.cesEst} returns a vector of the estimated coefficients.

   \code{coef.summary.cesEst} returns a matrix with four columns:
      the estimated coefficients/parameters of the CES,
      their standard errors, the t-statistic,
      and corresponding (two-sided) P-values.

   \code{vcov.cesEst} returns the variance covariance matrix
   of the estimated coefficients.
}

\seealso{\code{\link{cesEst}} and \code{\link{summary.cesEst}}.}

\author{Arne Henningsen and Geraldine Henningsen}

\examples{
   data( germanFarms )
   # output quantity:
   germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
   # quantity of intermediate inputs
   germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput


   ## CES: Land & Labor
   cesLandLabor <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms )

   # estimated coefficients
   coef( cesLandLabor )

   # estimated coefficients, their standard errors, t-statistic, P-values
   coef( summary( cesLandLabor ) )

   # covariance matrix of the estimated coefficients
   vcov( cesLandLabor )
}

\keyword{models}
