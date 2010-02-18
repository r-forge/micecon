\name{plot.cesEstGridRho}
\alias{plot.cesEstGridRho}
\encoding{UTF-8}
\title{Plot RSSs of a CES function Estimated by Grid Search}

\description{
   Plot a scatter plot,
   where the values of \eqn{\rho}{rho} are on the x axis
   and the corresponding sums of the squared residuals
   obtained by a grid search for \eqn{\rho}{rho}
   are on the y axis.
   Estimations that did not converge are marked with red.
}

\usage{
\method{plot}{cesEstGridRho}( x, \dots )
}

\arguments{
   \item{x}{object returned by \code{cesEst}
      if it was called with argument \code{rho} set a vector
      containing more than one value for \eqn{\rho}{rho}
      so that a grid search was performed.}
   \item{\dots}{All further arguments are passed
      to \code{\link[graphics]{plot.default}}.}
}

\seealso{
   \code{\link{cesEst}}.
}

\author{Arne Henningsen}

\examples{
   data( germanFarms )
   # output quantity:
   germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
   # quantity of intermediate inputs
   germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput

   ## CES: Land & Intermediate Inputs
   cesLandInt <- cesEst( yName = "qOutput",
      xNames = c( "land", "qVarInput" ), data = germanFarms,
      rho = seq( from = -0.6, to = 0.9, by = 0.3 ) )

   # plot the rhos against the sum of squared residuals
   plot( cesLandInt ) 
}

\keyword{models}
\keyword{regression}
\keyword{nonlinear}
