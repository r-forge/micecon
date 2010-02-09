\name{cesEstGridRho}
\alias{cesEstGridRho}
\encoding{UTF-8}
\title{Estimate a CES function by a grid search for rho}

\description{
   Estimate a Constant-Elasticities-of-Scale (CES) function
   using a one-dimensional grid search for \eqn{\rho}{rho}.
   The CES ist estimated by least-squares using \code{\link{cesEst}},
   where \eqn{\rho}{rho} is fixed at different values.
   Finally the estimation with the \eqn{\rho}{rho},
   that results in the smallest sum of squared residuals is chosen.
}

\usage{
cesEstGridRho(from = -0.8, to = 4, by = 0.2, returnAll = FALSE, \dots )
}

\arguments{
   \item{from}{first (minimal) value of \eqn{\rho}{rho}
      used in the grid search.}
   \item{to}{last (maximal) value of \eqn{\rho}{rho}
      used in the grid search.}
   \item{by}{increment of the values of \eqn{\rho}{rho}
      used in the grid search.}
   \item{returnAll}{logical. Return all objects returned
      by \code{\link{cesEst}} during the grid search?}
   \item{\dots}{All further arguments are passed to \code{\link{cesEst}};
      please note that they have to be named.}
}

\value{
   \code{cesEst} returns a list of class \code{cesEstGridRho}
   (and \code{cesEst})
   that is the object returned by \code{\link{cesEst}}
   with the value of \eqn{\rho}{rho}
   that results in the smallest sum of squared residuals.
   However, the returned object has one or two further components:
   \item{otherSum}{data frame with summary results of the estimations
      with the \sQuote{other} values of \eqn{\rho}{rho},
      i.e.\ the values that did \emph{not} result
      in the smalles sum of squared residuals;
      this data frame has follwing columns:
      \code{rho} = the value of \eqn{\rho}{rho},
      \code{rss} = the corresponding sum of squared residuals, and
      \code{convergence} = logical value indicating whether the estimation
      converged (if appropriate for the method used for the estimation).}
   \item{otherFull}{list of objects returned by \code{\link{cesEst}}
      with the \sQuote{other} values of \eqn{\rho}{rho},
      i.e.\ the values that did \emph{not} result
      in the smalles sum of squared residuals
      (only if argument \code{returnAll} is set to \code{TRUE}).}
}

\seealso{
   \code{\link{cesEst}}.
}

\author{Arne Henningsen and Geraldine Henningsen}

\examples{
   data( germanFarms )
   # output quantity:
   germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
   # quantity of intermediate inputs
   germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput

   ## CES: Land & Intermediate Inputs
   cesLandInt <- cesEstGridRho( yName = "qOutput",
      xNames = c( "land", "qVarInput" ), data = germanFarms,
      from = -0.6, to = 0.9, by = 0.3 )
}

\keyword{models}
\keyword{regression}
\keyword{nonlinear}
