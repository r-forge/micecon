\name{cesEst}
\alias{cesEst}
\alias{print.cesEst}
\title{Estimate a CES function}

\description{
   Estimate a Constant-Elasticities-of-Scale (CES) function
   with two exogenous variables by Least Squares.

   \deqn{y = \gamma * ( \delta * x1^{-\rho} + ( 1 - \delta ) * x2^{-\rho} )^{-\phi/\rho}}{%
      y = gamma * ( delta * x1^(-rho) + ( 1 - delta ) * x2^(-rho) )^(-phi/rho)}
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
cesEst( yName, xNames, data, vrs = FALSE, method = "Nelder-Mead",
   startVal = c( sqrt( mean( data[[ yName ]] ) ), 0.5, -0.5, 1 ),
   \dots )
\method{print}{cesEst}( x, digits = max(3, getOption("digits") - 3),
   \dots )
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
      estimation (using \code{\link[stats]{optim}}),
      or \code{"LM"} for non-linear least-squares using the Levenberg-Marquardt
      algorithm (using \code{\link[minpack.lm]{nls.lm}}).}
   \item{startVal}{starting values for the non-linear estimations:
      a numeric vector of length 3 for constant returns to scale
      (CRS, without \eqn{\phi}{phi}, i.e. \eqn{\phi}{phi} = 1)
      or of length 4 for variable returns to scale (VRS, with \eqn{\phi}{phi});
      the order must be \eqn{\gamma}{gamma}, \eqn{\delta}{delta},
      \eqn{\rho}{rho}, and eventually \eqn{\phi}{phi};
      names of the elements are ignored.}
   \item{x}{an object of class \code{cesEst}.}
   \item{digits}{number of digits.}
   \item{\dots}{further arguments to \code{cesEst} are passed to
      \code{\link[stats]{optim}} or \code{\link[minpack.lm]{nls.lm}};
      further arguments to \code{print.cesEst} are currently ignored.}
}

\details{
   Analytical gradients are used (only) if the \dQuote{BFGS}, \dQuote{CG},
   or \dQuote{L-BFGS-B} method is used.
}

\value{
   \code{cesEst} returns a list of class \code{cesEst}
   that contains following elements:
   \item{coefficients}{estimated coefficients/parameters of the CES.}
   \item{vcov}{approximate covariance matrix of the estimated parameters
      calculated from the parameters of the linearized model by the Delta method
      (only if argument \code{method} is \code{"Kmenta"}).}
   \item{cov.unscaled}{unscaled covariance matrix of the estimated parameters,
      i.e. the product of the gradient matrices evaluated at the estimated
      parameters.}
   \item{fitted.values}{the fitted values.}
   \item{residuals}{the residuals.}
   \item{call}{the matched call.}
   \item{method}{argument \code{method}.}
   \item{optim}{object returned by \code{\link[stats]{optim}}
      (only if argument \code{method} is \code{"Nelder-Mead"}, \code{"BFGS"},
      \code{"CG"}, \code{"L-BFGS-B"}, or \code{"SANN")}).}
   \item{nls.lm}{object returned by \code{\link[minpack.lm]{nls.lm}}
      (only if argument \code{method} is \code{"LM"}).}
   \item{translogUnr}{unrestricted estimation of the translog model
      (only if argument \code{method} is \code{"Kmenta"}).}
   \item{kmenta}{estimation results of the Kmenta approximation
      (a restricted translog model)
      returned by \code{\link[systemfit]{systemfit}}
      (only if argument \code{method} is \code{"Kmenta"}).}
   \item{test}{test for the restrictions of translog model
      (only if argument \code{method} is \code{"Kmenta"}).}
}

\seealso{
   \code{\link{summary.cesEst}} for the \code{summary} method,
   \code{\link{coef.cesEst}} for several further methods,
   \code{\link{cesCalc}} for calculations or simulations with the CES,
   \code{\link[micEcon]{translogEst}} for estimating translog functions, and
   \code{\link[micEcon]{quadFuncEst}} for estimating quadratic functions.}

\author{Arne Henningsen and Geraldine Henningsen}

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

   # using the L-BFGS-B optimization method with constrained delta
   cesLandLaborBfgsCon <- cesEst( "qOutput", c( "land", "qLabor" ),
      germanFarms, method = "L-BFGS-B", lower = c( -Inf, 0, -Inf ),
      upper = c( Inf, 1, Inf ) )

   # using the CG optimization method
   cesLandLaborSann <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
      method = "CG" )

   # using the SANN optimization method
   cesLandLaborSann <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
      method = "SANN" )

   # using the Kmenta approximation
   cesLandLaborSann <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
      method = "Kmenta" )


   ## CES: Land & Intermediate Inputs
   cesLandInt <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms )

   # variable returns to scale
   cesLandIntVrs <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
      vrs = TRUE )

   # using the BFGS optimization method
   cesLandIntBfgs <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
      method = "BFGS" )

   # using the L-BFGS-B optimization method with constrained delta
   cesLandIntBfgsCon <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
      method = "L-BFGS-B", lower = c( -Inf, 0, -Inf ),
      upper = c( Inf, 1, Inf ) )

   # using the SANN optimization method
   cesLandIntSann <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
      method = "SANN" )

   # using the Levenberg-Marquardt optimization method
   cesLandIntLm <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
      method = "LM", control = nls.lm.control( maxiter = 200 ) )
}

\keyword{models}
