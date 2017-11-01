\name{durbinWatsonTest.cesEst}
\alias{durbinWatsonTest.cesEst}
\alias{dwt.cesEst}
\encoding{UTF-8}
\title{Durbin Watson Test for Estimated CES Functions}

\description{
   Conduct a generalized Durbin-Watson-Test
   as suggested by White (1992)
   to test for serial correlation of the residuals.
   \code{dwt} is an abbreviation for \code{durbinWatsonTest}.
}

\usage{
\method{durbinWatsonTest}{cesEst}( model, \dots )
}

\arguments{
   \item{model}{object returned by \code{cesEst};
      in the call to \code{cesEst},
      argument \code{returnGrad} must be set to \code{TRUE}.}
   \item{\dots}{further arguments are passed
      to \code{\link[car]{durbinWatsonTest.lm}}.}
}

\references{
White, K.J. (1992):
The Durbin-Watson Test for Autocorrelation in Nonlinear Models.
The Review of Economics and Statistics 74(2), p. 370-373.
}

\seealso{
   \code{\link{cesEst}}, \code{\link[car]{durbinWatsonTest}}.
}

\author{Arne Henningsen}

\examples{
   data( germanFarms, package = "micEcon" )
   # output quantity:
   germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
   # quantity of intermediate inputs
   germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput

   ## CES: Land & Intermediate Inputs
   cesLandInt <- cesEst( yName = "qOutput",
      xNames = c( "land", "qVarInput" ), data = germanFarms,
      returnGrad = TRUE )

   # conduct the generalized Durbin-Watson test
   dwt( cesLandInt ) 
}

\keyword{models}
\keyword{regression}
\keyword{nonlinear}
