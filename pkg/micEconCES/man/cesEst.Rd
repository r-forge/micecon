\name{cesEst}
\alias{cesEst}
\title{Estimate a CES function}

\description{
   Estimate a Constant-Elasticities-of-Scale (CES) function
   with two exogenous variables by Least Squares.

   \deqn{y = \gamma * ( \alpha * x1^{-\rho} + ( 1 - \alpha ) * x2^{-\rho} )^{-\phi/\rho}}{%
      y = gamma * ( alpha * x1^(-rho) + ( 1 - alpha ) * x2^(-rho) )^(-phi/rho)}
   where the elesticity of substitution is
   \eqn{s = 1 / ( 1 + \rho )}{s = 1 / ( 1 + rho )}

   Warning: The econometric estimation of a CES function is (almost)
   always very problematic,
   because very different parameter vectors could result in very similar
   values of the objective function (sum of squared residuals).
   Hence, even if \code{\link[stats]{optim}} reports
   that the nonlinear minimization has converged,
   there might be another rather different parameter vector
   that results in a lower sum of squared residuals.
}

\usage{
cesEst( yName, xNames, data, vrs = FALSE, method = "Nelder-Mead", \dots )
}

\arguments{
   \item{yName}{a string containing the name of the dependent
      variable.}
   \item{xNames}{a vector of two strings containing the names of the
      independent variables.}
   \item{data}{data frame containing the data.}
   \item{vrs}{logical. Allow for variable returns to scale?}
   \item{method}{character string indicationg the estimation method:
      either \code{"Kmenta"} for the Kmenta approximation
      or \code{"Nelder-Mead"}, \code{"BFGS"}, \code{"CG"},
      \code{"L-BFGS-B"}, or \code{"SANN")} for a non-linear least-squares
      estimation (using \code{\link[stats]{optim}}).}
   \item{\dots}{further arguments are passed to \code{\link[stats]{optim}}.}
}

\details{
   Analytical gradients are used (only) if the \dQuote{BFGS}, \dQuote{CG},
   or \dQuote{L-BFGS-B} method is used.
}

\value{
   A list of class \code{cesEst}
   that contains following elements:
   \item{coefficients}{estimated coefficients/parameters of the CES.}
   \item{fitted.values}{the fitted values.}
   \item{vcov}{covariance matrix of the estimated parameters.}
   \item{call}{the matched call.}
   \item{method}{argument \code{method}.}
   \item{optim}{object returned by \code{\link[stats]{optim}}
      (only if argument \code{method} is \code{"Nelder-Mead"}, \code{"BFGS"},
      \code{"CG"}, \code{"L-BFGS-B"}, or \code{"SANN")}).}
   \item{translogUnr}{unrestricted estimation of the translog model
      (only if argument \code{method} is \code{"Kmenta"}).}
   \item{translog}{restricted estimation of the translog model
      (only if argument \code{method} is \code{"Kmenta"}).}
   \item{test}{test for the restrictions of translog model
      (only if argument \code{method} is \code{"Kmenta"}).}
}

\seealso{\code{\link[micEcon]{translogEst}} and \code{\link[micEcon]{quadFuncEst}}.}

\author{Arne Henningsen}

\examples{
   data( germanFarms )
   # output quantity:
   germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
   # quantity of intermediate inputs
   germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput


   ## CES: Land & Labor
   cesLandLabor <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms )

   # variable returns to scale
   cesLandLaborVrs <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
      vrs = TRUE )

   # using the BFGS optimization method
   cesLandLaborBfgs <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
      method = "BFGS" )

   # using the L-BFGS-B optimization method with constrained alpha
   cesLandLaborBfgsCon <- cesEst( "qOutput", c( "land", "qLabor" ),
      germanFarms, method = "L-BFGS-B", lower = c( -Inf, 0, -Inf ),
      upper = c( Inf, 1, Inf ) )

   # using the SANN optimization method
   cesLandLaborSann <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
      method = "SANN" )


   ## CES: Land & Intermediate Inputs
   cesLandInt <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms )

   # variable returns to scale
   cesLandIntVrs <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
      vrs = TRUE )

   # using the BFGS optimization method
   cesLandIntBfgs <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
      method = "BFGS" )

   # using the L-BFGS-B optimization method with constrained alpha
   cesLandIntBfgsCon <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
      method = "L-BFGS-B", lower = c( -Inf, 0, -Inf ),
      upper = c( Inf, 1, Inf ) )

   # using the SANN optimization method
   cesLandIntSann <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
      method = "SANN" )
}

\keyword{models}