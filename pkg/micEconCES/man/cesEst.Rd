\name{cesEst}
\alias{cesEst}
\alias{print.cesEst}
\encoding{UTF-8}
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
   Hence, even if the optimizer reports
   that the nonlinear minimization has converged,
   there might be another rather different parameter vector
   that results in a lower sum of squared residuals.
}

\usage{
cesEst( yName, xNames, data, vrs = FALSE, method = "Nelder-Mead",
   startVal = NULL, lower = "auto", upper = "auto", \dots )
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
      \code{"L-BFGS-B"}, \code{"SANN"}, \code{"LM"},
      \code{"Newton"}, \code{"PORT"}, or \code{"DE"}
      for non-linear least-squares (see section \sQuote{Details}).}
   \item{startVal}{optional numeric vector giving the starting values
      of the parameters in the non-linear estimations
      (see section \sQuote{Details}).}
   \item{lower}{lower bounds of the parameters
      (see section \sQuote{Details}).}
   \item{upper}{upper bounds of the parameters
      (see section \sQuote{Details}).}
   \item{x}{an object of class \code{cesEst}.}
   \item{digits}{number of digits.}
   \item{\dots}{further arguments to \code{cesEst} are passed to
      \code{\link[stats]{optim}}, \code{\link[minpack.lm]{nls.lm}},
      \code{\link[stats]{nlm}}, \code{\link[stats]{nlminb}},
      or \code{\link[DEoptim]{DEoptim}};
      further arguments to \code{print.cesEst} are currently ignored.}
}

\details{
   \bold{Estimation method}\cr
   Argument \code{method} determines the estimation method.
   If it is \code{"Kmenta"},
   the CES is estimated by ordinary least squares using the Kmenta
   approximation;
   otherwise, it is estimated by non-linear least-squares.
   Several different optimizers can be used for the non-linear estimation.
   The optimization methods
   \code{Nelder-Mead} (Nelder and Mead 1965),
   \code{BFGS} (Broyden 1970, Fletcher 1970, Goldfarb 1970, Shanno 1970),
   \code{CG} (Conjugate Gradients based on Fletcher and Reeves 1964),
   \code{L-BFGS-B} (with box-constraints, Byrd, Lu, Nocedal, and Zhu 1995), and
   \code{SANN} (Simulated Annealing, Bélisle 1992)
   use \code{\link[stats]{optim}} for the optimization.
   The optimization method
   \code{LM} (Levenberg-Marquardt, see Moré 1978)
   uses \code{\link[minpack.lm]{nls.lm}} for the optimization.
   The optimization method
   \code{Newton} (Newton-type, see Dennis and Schnabel 1983
   and Schnabel, Koontz, and Weiss 1985)
   uses \code{\link[stats]{nlm}} for the optimization.
   The optimization method
   \code{PORT} (PORT routines, see Gay 1990)
   uses \code{\link[stats]{nlminb}} for the optimization.
   The optimization method
   \code{DE} (Differential Evolution, see Storn and Price 1997)
   uses \code{\link[DEoptim]{DEoptim}} for the optimization.
   Analytical gradients are used in the \code{BFGS}, \code{CG},
   \code{L-BFGS-B}, \code{LM}, \code{Newton}, and \code{PORT} method.\cr\cr


   \bold{Starting values}\cr
   Argument \code{startVal} should be a numeric vector
   of length 3 for estimations of the CES with constant returns to scale
   (CRS, i.e. with \eqn{\phi}{phi} fixed at 1)
   or of length 4 for estimations of the CES with variable returns to scale
   (VRS, i.e. with \eqn{\phi}{phi} estimated).
   The order must be \eqn{\gamma}{gamma}, \eqn{\delta}{delta},
   \eqn{\rho}{rho}, and eventually \eqn{\phi}{phi}.
   Names of the elements are ignored.
   If argument \code{startVal} is \code{NULL},
   the starting value of \eqn{\delta}{delta} is set to 0.5,
   the starting value of \eqn{\rho}{rho} is set to 0.25
      (i.e.\ elasticity of substitution = 0.8),
   the starting value of \eqn{\phi}{phi} is set to 1
      (only if argument \code{vrs} is \code{TRUE}),
   and \eqn{\gamma}{gamma} is set to a value
   so that the mean of the endogenous variable
   is equal to the mean of its fitted values,
   i.e. \eqn{\gamma}{gamma} =  mean( y ) /
      mean( CES( X, \eqn{\gamma}{gamma} = 1, \eqn{\delta}{delta} = 0.5,
         \eqn{\rho}{rho} = 0.25, \eqn{\phi}{phi} = 1 ) ).\cr\cr

   \bold{Lower and upper bounds}\cr
   Arguments \code{lower} and \code{upper} can be used
   to set lower and upper bounds for the estimated parameters.
   If these arguments are \code{-Inf} and \code{Inf}, respectively,
   the parameters are estimated without unconstraints.
   By default, arguments \code{lower} and \code{upper} are both equal
   to the character string \code{"auto"},
   which means that the bounds are set automatically
   depending on the estimation method:
   In case of the \code{L-BFGS-B}, \code{PORT}, and \code{DE} method,
   the lower bound is \code{0} for \eqn{\gamma}{gamma},
   \code{0} for \eqn{\delta}{delta}, \code{-1} for \eqn{\rho}{rho},
   and eventually \code{0} for \eqn{\phi}{phi}.
   In case of the \code{L-BFGS-B} and \code{PORT} method,
   the upper bound is infinity for \eqn{\gamma}{gamma},
   \code{1} for \eqn{\delta}{delta}, infinity for \eqn{\rho}{rho},
   and eventually infinity for \eqn{\phi}{phi}.
   Since the \sQuote{Differential Evulation} algorithm requires finit bounds,
   the upper bounds for the \code{DE} method are set
   to \code{1e10} for \eqn{\gamma}{gamma},
   \code{1} for \eqn{\delta}{delta}, \code{10} for \eqn{\rho}{rho},
   and eventually \code{10} for \eqn{\phi}{phi}.
   In case of all other estimation methods,
   the lower and upper bounds are set to \code{-Inf} and \code{Inf}, respectively,
   because these methods do not support parameter constraints.
   Of course, the user can specify own lower and upper bounds
   by setting arguments \code{lower} and \code{upper} to numeric vectors
   that should have the same format as argument \code{startVal}
   (see above).
}

\value{
   \code{cesEst} returns a list of class \code{cesEst}
   that has following components:
   \item{coefficients}{estimated coefficients/parameters of the CES.}
   \item{iter}{number of iterations
      (only for non-linear least-squares estimations).}
   \item{convergence}{logical value indicating if the non-linear estimation
      has converged (only for non-linear least-squares estimations).}
   \item{message}{additional information from the optimizer
      (only if a message was returned by \code{\link{optim}} or
      \code{\link[minpack.lm]{nls.lm}}.}
   \item{vcov}{approximate covariance matrix of the estimated parameters
      calculated from the parameters of the linearized model by the Delta method
      (only if argument \code{method} is \code{"Kmenta"}).}
   \item{cov.unscaled}{unscaled covariance matrix of the estimated parameters,
      i.e. the inverse of the cross-product of the gradient matrix
      evaluated at the estimated parameters.}
   \item{fitted.values}{the fitted values.}
   \item{residuals}{the residuals.}
   \item{call}{the matched call.}
   \item{method}{argument \code{method}.}
   \item{startVal}{starting values for the non-linear estimation
      (not for the \code{Kmenta} and \code{DE} method).}
   \item{lower}{lower bounds of the parameters.}
   \item{upper}{upper bounds of the parameters.}
   \item{optim}{object returned by \code{\link[stats]{optim}}
      (only if argument \code{method} is \code{"Nelder-Mead"}, \code{"BFGS"},
      \code{"CG"}, \code{"L-BFGS-B"}, or \code{"SANN"})).}
   \item{nls.lm}{object returned by \code{\link[minpack.lm]{nls.lm}}
      (only if argument \code{method} is \code{"LM"}).}
   \item{nlm}{object returned by \code{\link[stats]{nlm}}
      (only if argument \code{method} is \code{"Newton"}).}
   \item{nlminb}{object returned by \code{\link[stats]{nlminb}}
      (only if argument \code{method} is \code{"PORT"}).}
   \item{DEoptim}{object returned by \code{\link[DEoptim]{DEoptim}}
      (only if argument \code{method} is \code{"DE"}).}
   \item{translog}{estimation results of the (unrestricted) translog model
      returned by \code{\link[micEcon]{translogEst}}
      (only if argument \code{method} is \code{"Kmenta"}).}
   \item{kmenta}{estimation results of the Kmenta approximation
      (a restricted translog model)
      returned by \code{\link[systemfit]{systemfit}}
      (only if argument \code{method} is \code{"Kmenta"}).}
   \item{testKmenta}{test of the restrictions implied
      by the Kmenta approximation in the unrestricted translog model
      returned by \code{\link[car]{linear.hypothesis}}
      (only if argument \code{method} is \code{"Kmenta"}).}
   \item{testCrs}{test of the restrictions of constant returns to scale
      (these restrictions include the restrictions implied
      by the Kmenta approximation)
      in the unrestricted translog model
      returned by \code{\link[car]{linear.hypothesis}}
      (only if argument \code{vrs} is \code{FALSE}
      and argument \code{method} is \code{"Kmenta"}).}
}

\references{
   Bélisle, C.J.P. (1992):
   Convergence theorems for a class of simulated annealing algorithms on Rd,
   Journal of Applied Probability 29, p. 885-895.

   Broyden, C.G. (1970):
   The Convergence of a Class of Double-rank Minimization Algorithms,
   Journal of the Institute of Mathematics and Its Applications 6, p. 76-90.

   Byrd, R.H., Lu, P., Nocedal, J. and Zhu, C. (1995):
   A limited memory algorithm for bound constrained optimization,
   SIAM J. Scientific Computing 16, p. 1190-1208.

   Dennis, J.E. and Schnabel, R.B. (1983):
   Numerical Methods for Unconstrained Optimization and Nonlinear Equations,
   Prentice-Hall, Englewood Cliffs, NJ.

   Fletcher, R. (1970):
   A New Approach to Variable Metric Algorithms,
   Computer Journal 13, p. 317-322.

   Fletcher, R. and Reeves, C.M. (1964):
   Function minimization by conjugate gradients,
   Computer Journal 7, p. 148-154.

   Gay, D.M. (1990):
   Usage Summary for Selected Optimization Routines,
   Computing Science Technical Report No. 153,
   AT&T Bell Laboratories, Murray Hill NJ,
   \url{http://netlib.bell-labs.com/cm/cs/cstr/153.pdf}.

   Goldfarb, D. (1970):
   A Family of Variable Metric Updates Derived by Variational Means,
   Mathematics of Computation 24, p. 23-26.

   Moré, J.J. (1978):
   The Levenberg-Marquardt algorithm: implementation and theory,
   in G.A. Watson (Ed.), Lecture Notes in Mathematics 630: Numerical Analysis,
   pp. 105-116, Springer-Verlag: Berlin.

   Nelder, J.A. and Mead, R. (1965):
   A simplex algorithm for function minimization,
   Computer Journal 7, p. 308-313.

   Schnabel, R.B., Koontz, J.E. and Weiss, B.E. (1985):
   A modular system of algorithms for unconstrained minimization,
   ACM Trans. Math. Software, 11, pp. 419-440.

   Shanno, D.F. (1970):
   Conditioning of Quasi-Newton Methods for Function Minimization,
   Mathematics of Computation 24, p. 647-656.

   Storn, R. and Price, K. (1997):
   Differential Evolution - A Simple and Efficient Heuristic
   for Global Optimization over Continuous Spaces,
   Journal of Global Optimization, 11(4), p. 341-359.
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

   # using the L-BFGS-B optimization method with constrained parameters
   cesLandLaborBfgsCon <- cesEst( "qOutput", c( "land", "qLabor" ),
      germanFarms, method = "L-BFGS-B" )

   # using the CG optimization method
   cesLandLaborSann <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
      method = "CG" )

   # using the SANN optimization method
   cesLandLaborSann <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
      method = "SANN" )

   # using the Kmenta approximation
   cesLandLaborKmenta <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
      method = "Kmenta" )

   # using the PORT optimization routine with unconstrained parameters
   cesLandLaborPortCon <- cesEst( "qOutput", c( "land", "qLabor" ),
      germanFarms, vrs = TRUE, method = "PORT", lower = -Inf, upper = Inf )

   # using the PORT optimization routine with constrained parameters and VRS
   cesLandLaborPortCon <- cesEst( "qOutput", c( "land", "qLabor" ),
      germanFarms, vrs = TRUE, method = "PORT" )

   # using the Differential Evolution optimization method
   cesLandLaborDe <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
      method = "DE" )


   ## CES: Land & Intermediate Inputs
   cesLandInt <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms )

   # variable returns to scale
   cesLandIntVrs <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
      vrs = TRUE )

   # using the BFGS optimization method
   cesLandIntBfgs <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
      method = "BFGS" )

   # using the L-BFGS-B optimization method with constrained parameters, VRS
   cesLandIntBfgsCon <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
      method = "L-BFGS-B", vrs = TRUE )

   # using the SANN optimization method
   cesLandIntSann <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
      method = "SANN" )

   # using the Levenberg-Marquardt optimization method
   cesLandIntLm <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
      method = "LM", control = nls.lm.control( maxiter = 200 ) )

   # using the Differential Evolution optimization method
   cesLandIntDe <- cesEst( "qOutput", c( "land", "qVarInput" ), germanFarms,
      method = "DE" )
}

\keyword{models}
\keyword{regression}
\keyword{nonlinear}
