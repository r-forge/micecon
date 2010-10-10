\name{cesEst}
\alias{cesEst}
\alias{print.cesEst}
\encoding{UTF-8}
\title{Estimate a CES function}

\description{
   Estimate a Constant-Elasticities-of-Scale (CES) function
   with two exogenous variables
   or a nested Constant-Elasticities-of-Scale (CES) function
   proposed by Sato (1967) with three or four exogenous variables
   by Least Squares.
   The functional forms are shown in the documentation of
   function \code{\link{cesCalc}}.

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
cesEst( yName, xNames, data, vrs = FALSE, method = "LM",
   start = NULL, lower = NULL, upper = NULL, rho = NULL,
   returnGridAll = FALSE, random.seed = 123,
   rhoApprox = c( y = 5e-6, gamma = 5e-6, delta = 5e-6, 
      rho = 1e-3, nu = 5e-6 ),
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
      or \code{"LM"}, \code{"NM"}, \code{"Nelder-Mead"}, \code{"BFGS"},
      \code{"CG"}, \code{"L-BFGS-B"}, \code{"SANN"},
      \code{"Newton"}, \code{"PORT"}, or \code{"DE"}
      for non-linear least-squares (see section \sQuote{Details}).}
   \item{start}{optional numeric vector giving the starting values
      of the parameters in the non-linear estimations
      (see section \sQuote{Details}).}
   \item{lower}{lower bounds of the parameters
      (see section \sQuote{Details}).}
   \item{upper}{upper bounds of the parameters
      (see section \sQuote{Details}).}
   \item{rho}{numeric scalar or vector at which the coefficient \eqn{\rho}{rho}
      should be fixed;
      if \code{NULL} (default), \eqn{\rho}{rho} is estimated together
      with the other parameters;
      if it has more than one element,
      a one-dimensional grid search for \eqn{\rho}{rho} is performed
      (see section \sQuote{Details}).}
   \item{returnGridAll}{logical value that indicates
      whether the estimates for all values of \eqn{\rho}{rho}
      obtained during the grid search
      (not just the estimations with the \sQuote{best} \eqn{\rho}{rho})
      should be returned
      (ignored if argument \code{rho} is \code{NULL}
      or has only a single element).}
   \item{random.seed}{an integer used to seed R's random number generator.
      This is to ensure replicability when the \code{"SANN"}
      or \code{"DE"} method is used. Defaults to 123.}
   \item{rhoApprox}{numeric vector with exactly 5 elements;
      the endogenous variable of the CES
      and the derivatives with respect to its coefficients
      are calculated using a first-order Taylor series approximation
      at \eqn{\rho=0}{rho = 0} (non-nested CES)
      or by interpolation between \eqn{\rho}{rho}, \eqn{\rho_1}{rho_1}, 
      or \eqn{\rho_2}{rho_2} equal to zero 
      and \eqn{\rho}{rho}, \eqn{\rho_1}{rho_1}, or \eqn{\rho_2}{rho_2} 
      equal to \eqn{\pm}{+/-}\code{rhoApprox}
      (nested CES),
      if the absolute value of the coefficients \eqn{\rho}{rho}, 
      \eqn{\rho_1}{rho_1}, or \eqn{\rho_2}{rho_2}
      is smaller than or equal to the corresponding element of this argument
      (see also argument \code{rhoApprox} of \code{\link{cesCalc}});
      the first element determines the threshold for calculating
      the endogenous variable; 
      the second element determines the threshold for calculating
      the derivatives with respect to \eqn{\gamma}{gamma},
      \eqn{\gamma_1}{gamma_1}, and \eqn{\gamma_2}{gamma_2};
      the third element determines the threshold for calculating
      the derivatives with respect to \eqn{\delta}{delta},
      \eqn{\delta_1}{delta_1}, and \eqn{\delta_2}{delta_2};
      the fourth element determines the threshold for calculating
      the derivatives with respect to \eqn{\rho}{rho},
      \eqn{\rho_1}{rho_1}, and \eqn{\rho_2}{rho_2};
      the fifth element determines the threshold for calculating
      the derivatives with respect to \eqn{\nu}{nu}.}
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
   The optimization method
   \code{LM} (Levenberg-Marquardt, see Moré 1978)
   uses \code{\link[minpack.lm]{nls.lm}} for the optimization.
   The optimization methods
   \code{NM} or \code{Nelder-Mead} (Nelder and Mead 1965),
   \code{BFGS} (Broyden 1970, Fletcher 1970, Goldfarb 1970, Shanno 1970),
   \code{CG} (Conjugate Gradients based on Fletcher and Reeves 1964),
   \code{L-BFGS-B} (with box-constraints, Byrd, Lu, Nocedal, and Zhu 1995), and
   \code{SANN} (Simulated Annealing, Bélisle 1992)
   use \code{\link[stats]{optim}} for the optimization.
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
   Analytical gradients are used in the \code{LM}, \code{BFGS}, \code{CG},
   \code{L-BFGS-B}, \code{Newton}, and \code{PORT} method.\cr


   \bold{Starting values}\cr
   Argument \code{start} should be a numeric vector.
   The order must be as described in the documentation of argument \code{coef}
   of function \code{\link{cesCalc}}.
   However, names of the elements are ignored.
   If argument \code{start} is \code{NULL},
   pre-defined starting values are used.
   The starting values of \eqn{\gamma_1}{gamma_1} (if present) is set to 1,
   the starting values of \eqn{\delta}{delta}, \eqn{\delta_1}{delta_1}, 
   and \eqn{\delta_2}{delta_2} (if present) are set to 0.5,
   the starting values of \eqn{\rho_1}{rho_1}, \eqn{\rho_2}{rho_2}, and
   \eqn{\rho}{rho} (if present and required) are set to 0.25
      (i.e.\ elasticity of substitution = 0.8 in the two-input case),
   the starting value of \eqn{\nu}{nu} (if required) is set to 1,
   and the starting value of \eqn{\gamma}{gamma} or \eqn{\gamma_2}{gamma_2} 
   is set to a value
   so that the mean of the endogenous variable
   is equal to the mean of its fitted values.
   Hence, for the CES with two inputs
   \eqn{\gamma}{gamma} =  mean( y ) /
      mean( CES( X, \eqn{\gamma}{gamma} = 1, \eqn{\delta}{delta} = 0.5,
         \eqn{\rho}{rho} = fixed or 0.25, \eqn{\nu}{nu} = 1 ) ),
   for the nested CES with 3 inputs
   \eqn{\gamma_2}{gamma_2} =  mean( y ) /
      mean( CES( X, \eqn{\gamma_1}{gamma_1} = 1, \eqn{\gamma_2}{gamma_2} = 1,
         \eqn{\delta_1}{delta_1} = 0.5, \eqn{\delta_2}{delta_2} = 0.5, 
         \eqn{\rho_1}{rho_1} = 0.25, \eqn{\rho}{rho} = fixed or 0.25, 
          \eqn{\nu}{nu} = 1 ) ),
   and for the nested CES with 4 inputs
   \eqn{\gamma}{gamma} =  mean( y ) /
      mean( CES( X, \eqn{\gamma}{gamma} = 1, \eqn{\delta_1}{delta_1} = 0.5,
          \eqn{\delta_2}{delta_2} = 0.5, \eqn{\rho_1}{rho_1} = 0.25, 
          \eqn{\rho_2}{rho_2} = 0.25, \eqn{\rho}{rho} = fixed or 0.25, 
          \eqn{\nu}{nu} = 1 ) ).\cr

   \bold{Lower and upper bounds}\cr
   Arguments \code{lower} and \code{upper} can be used
   to set lower and upper bounds for the estimated parameters.
   If these arguments are \code{-Inf} and \code{Inf}, respectively,
   the parameters are estimated without unconstraints.
   By default, arguments \code{lower} and \code{upper} are both \code{NULL},
   which means that the bounds are set automatically
   depending on the estimation method:
   In case of the \code{L-BFGS-B}, \code{PORT}, and \code{DE} method,
   the lower bound is \code{0} for \eqn{\gamma}{gamma}, \eqn{\gamma_1}{gamma_1}, 
   \eqn{\gamma_2}{gamma_2}, \eqn{\delta}{delta}, \eqn{\delta_1}{delta_1}, 
   and \eqn{\delta_2}{delta_2} (if present), 
   \code{-1} for \eqn{\rho_1}{rho_1}, \eqn{\rho_2}{rho_2}, and \eqn{\rho}{rho}
   (if present),
   and eventually \code{0} for \eqn{\nu}{nu}.
   In case of the \code{L-BFGS-B} and \code{PORT} method,
   the upper bound is infinity for \eqn{\gamma}{gamma}, \eqn{\gamma_1}{gamma_1}, 
   and \eqn{\gamma_2}{gamma_2} (if present),
   \code{1} for \eqn{\delta}{delta}, \eqn{\delta_1}{delta_1}, and
   \eqn{\delta_2}{delta_2} (if present),
   infinity for \eqn{\rho_1}{rho_1}, \eqn{\rho_2}{rho_2}, and \eqn{\rho}{rho}
   (if present),
   and eventually infinity for \eqn{\nu}{nu}.
   Since the \sQuote{Differential Evulation} algorithm requires finit bounds,
   the upper bounds for the \code{DE} method are set
   to \code{1e10} for \eqn{\gamma}{gamma}, \eqn{\gamma_1}{gamma_1},
   and \eqn{\gamma_2}{gamma_2} (if present),
   \code{1} for \eqn{\delta}{delta}, \eqn{\delta_1}{delta_1}, and 
   \eqn{\delta_2}{delta_2} (if present),
   \code{10} for \eqn{\rho_1}{rho_1}, \eqn{\rho_2}{rho_2}, and \eqn{\rho}{rho}
   (if present),
   and eventually \code{10} for \eqn{\nu}{nu}.
   In case of all other estimation methods,
   the lower and upper bounds are set to \code{-Inf} and \code{Inf}, respectively,
   because these methods do not support parameter constraints.
   Of course, the user can specify own lower and upper bounds
   by setting arguments \code{lower} and \code{upper} to numeric vectors
   that should have the same format as argument \code{start}
   (see above).\cr

   \bold{Grid search for \eqn{\rho}{rho}}\cr
   If argument \code{rho} has more than one element,
   a one-dimensional grid search for \eqn{\rho}{rho} is performed.
   The CES ist estimated by least-squares,
   where \eqn{\rho}{rho} is fixed consecutively at each value
   defined in argument \code{rho}.
   Finally the estimation with the \eqn{\rho}{rho},
   that results in the smallest sum of squared residuals is chosen
   (and returned).\cr

   \bold{Random numbers}\cr
   The \sQuote{state} (or \sQuote{seed}) of R's random number generator
   is saved at the beginning of the \code{cesEst} function
   and restored at the end of this function
   so that this function does \emph{not} affect the generation of random numbers
   although the random seed is set to argument \code{random.seed}
   and the \sQuote{SANN} and \sQuote{DE} algorithms use random numbers.
}

\value{
   \code{cesEst} returns a list of class \code{cesEst}
   that has following components:
   \item{coefficients}{estimated coefficients/parameters of the CES
      (including a possible fixed \eqn{\rho}{rho}).}
   \item{iter}{number of iterations
      (only for non-linear least-squares estimations).}
   \item{convergence}{logical value indicating if the non-linear estimation
      has converged (only for non-linear least-squares estimations
      with solvers that have a convergence criterion).}
   \item{message}{additional information from the optimizer
      (only if a message was returned by \code{\link{optim}} or
      \code{\link[minpack.lm]{nls.lm}}.}
   \item{vcov}{approximate covariance matrix of the estimated parameters
      calculated from the parameters of the linearized model by the Delta method
      (only if argument \code{method} is \code{"Kmenta"}).}
   \item{cov.unscaled}{unscaled covariance matrix of the estimated parameters
      (including a possible fixed \eqn{\rho}{rho}),
      i.e. the inverse of the cross-product of the gradient matrix
      evaluated at the estimated parameters.}
   \item{fitted.values}{the fitted values.}
   \item{residuals}{the residuals.}
   \item{rss}{the sum of the squared residuals 
      (i.e. the value of the objective function of the non-linear 
      least-squares estimation evaluated at the estimated parameters).}
   \item{call}{the matched call.}
   \item{method}{argument \code{method}.}
   \item{start}{starting values for the non-linear estimation
      (not for the \code{Kmenta} and \code{DE} method).}
   \item{lower}{lower bounds of the parameters.}
   \item{upper}{upper bounds of the parameters.}
   \item{rho}{argument \code{rho}.}
   \item{nls.lm}{object returned by \code{\link[minpack.lm]{nls.lm}}
      (only if argument \code{method} is \code{"LM"}).}
   \item{optim}{object returned by \code{\link[stats]{optim}}
      (only if argument \code{method} is \code{"NM"}, \code{"Nelder-Mead"},
      \code{"BFGS"}, \code{"CG"}, \code{"L-BFGS-B"}, or \code{"SANN"})).}
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
      by the Kmenta approximation
      (including constant returns to scale
      if argument \code{vrs} is \code{FALSE})
      in the unrestricted translog model
      returned by \code{\link[car]{linear.hypothesis}}
      (only if argument \code{method} is \code{"Kmenta"}).}
   \item{allRhoSum}{data frame with summary results of the estimations
      with all values of \eqn{\rho}{rho} used in the grid search
      (only if a grid search was performed);
      this data frame has follwing columns:
      \code{rho} = the value of \eqn{\rho}{rho},
      \code{rss} = the corresponding sum of squared residuals, and
      (if appropriate for the method used for the estimation)
      \code{convergence} = logical value indicating whether the estimation
      converged.}
   \item{allRhoFull}{list of estimation results returned by \code{\link{cesEst}}
      for all values of \eqn{\rho}{rho} used in the grid search
      (only if a grid search was performed and
      argument \code{returnGridAll} is set to \code{TRUE}).}
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
   \code{\link{plot.cesEst}} for plotting the results
      of the grid search for \eqn{\rho}{rho},
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


   ## CES: Land & Labor (Levenberg-Marquardt algorithm)
   cesLandLabor <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms )

   # variable returns to scale, increased max. number of iter. (LM algorithm)
   cesLandLaborVrs <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
      vrs = TRUE, control = nls.lm.control( maxiter = 1000 ) )

   # using the Nelder-Mead optimization method
   cesLandLaborNm <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
      method = "NM" )

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
   # (with decreased number of iteration to decrease execution time)
   cesLandLaborSann <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
      method = "SANN", control = list( maxit = 1000 ) )

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
   # (with decreased number of iteration to decrease execution time)
   cesLandLaborDe <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
      method = "DE", control = DEoptim.control( itermax = 50 ) )

   ## estimation with a grid search for rho (using the LM algorithm)
   cesLandInt <- cesEst( "qOutput", c( "land", "qLabor" ),
      data = germanFarms, rho = seq( from = -0.6, to = 0.9, by = 0.3 ) )
}

\keyword{models}
\keyword{regression}
\keyword{nonlinear}
