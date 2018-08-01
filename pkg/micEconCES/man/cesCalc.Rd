\name{cesCalc}
\alias{cesCalc}
\title{Calculate CES function}

\description{
   Calculate the endogenous variable
   of a \sQuote{Constant Elasticity of Substitution} (CES)
   function.

   The original CES function with two explanatory variables is

   \deqn{y = \gamma \: \exp ( \lambda \: t ) \: 
      ( \delta \: x_1^{-\rho} + ( 1 - \delta ) \: x_2^{-\rho} )
      ^{-\frac{\nu}{\rho}}}{%
      y = gamma * exp( lambda * t ) * 
         ( delta * x1^(-rho) + ( 1 - delta ) * x2^(-rho) )^(-nu/rho)}

   and the non-nested CES function with \eqn{N} explanatory variables is

   \deqn{y = \gamma \: \exp ( \lambda \: t ) \: 
      \left( \sum_{i=1}^N \delta_i \: x_i^{-\rho} \right)
      ^{-\frac{\nu}{\rho}}}{%
      y = gamma * exp( lambda * t ) * 
         ( sum(i=1 to N) delta_i * x_i^(-rho) )^(-nu/rho)}

   where in the latter case
   \eqn{\sum_{i=1}^N \delta_i = 1}{sum(i=1 to N) delta_i = 1}.

   In both cases, the elesticity of substitution is
   \eqn{s = \frac{1}{ 1 + \rho }}{s = 1 / ( 1 + rho )}.

   The \emph{nested} CES function with 3 explanatory variables 
   proposed by Sato (1967) is

   \deqn{y = \gamma \: \exp ( \lambda \: t ) \: \left[
      \delta \:
      \left( \delta_1 \: x_1^{-\rho_1} +
      ( 1 - \delta_1 ) x_2^{-\rho_1} \right)^{\frac{\rho}{\rho_1}} +
      ( 1 - \delta ) x_3^{-\rho}
      \right]^{-\frac{\nu}{\rho}}}{%
      y = gamma * exp( lambda * t ) * [ delta *
      ( delta_1 * x_1^(-rho_1) + ( 1 - delta_1 ) * x_2^(-rho_1) )^(rho / rho_1) +
      ( 1 - delta ) * x_3^(-rho)
      ]^(-nu / rho)}

   and the \emph{nested} CES function with 4 explanatory variables 
   (a generalisation of the version proposed by Sato, 1967) is

   \deqn{y = \gamma \: \exp ( \lambda \: t ) \: \left[ \delta \cdot
      \left( \delta_1 \: x_1^{-\rho_1} +
      ( 1 - \delta_1 ) x_2^{-\rho_1} \right)^{\frac{\rho}{\rho_1}} +
      ( 1 - \delta ) \cdot
      \left( \delta_2 \: x_3^{-\rho_2} +
      ( 1 - \delta_2 ) x_4^{-\rho_2} \right)^{\frac{\rho}{\rho_2}}
      \right]^{-\frac{\nu}{\rho}}}{%
      y = gamma * exp( lambda * t ) * 
      [ delta ( delta_1 * x_1^(-rho_1) + ( 1 - delta_1 ) * x_2^(-rho_1) )^(rho / rho_1) +
      ( 1 - delta ) ( delta_2 * x_3^(-\rho_2) + ( 1 - delta_2 ) * x_4^(-rho_2) )^(rho / rho_2)
      ]^(-nu / rho)}
}

\usage{
cesCalc( xNames, data, coef, tName = NULL, nested = FALSE, rhoApprox = 5e-6 )
}

\arguments{
   \item{xNames}{a vector of strings containing the names of the
      explanatory variables.}
   \item{data}{data frame containing the explanatory variables.}
   \item{coef}{numeric vector containing the coefficients of the CES:
      if the vector is unnamed,
      the order of the coefficients must be
      \eqn{\gamma}{gamma}, eventuelly \eqn{\lambda}{lambda}, 
      \eqn{\delta}{delta}, \eqn{\rho}{rho}, 
      and eventually \eqn{\nu}{nu}
      in case of two expanatory variables,
      \eqn{\gamma}{gamma}, eventuelly \eqn{\lambda}{lambda}, 
      \eqn{\delta_1}{delta_1}, \ldots,
      \eqn{\delta_N}{delta_N}, \eqn{\rho}{rho}, and eventually \eqn{\nu}{nu}
      in case of the non-nested CES with \eqn{N>2} explanatory variables,
      \eqn{\gamma}{gamma}, eventuelly \eqn{\lambda}{lambda}, 
      \eqn{\delta_1}{delta_1}, \eqn{\delta}{delta}, 
      \eqn{\rho_1}{rho_1}, \eqn{\rho}{rho}, and eventually \eqn{\nu}{nu}
      in case of the nested CES with 3 explanatory variables,
      and \eqn{\gamma}{gamma}, eventuelly \eqn{\lambda}{lambda}, 
      \eqn{\delta_1}{delta_1}, \eqn{\delta_2}{delta_2}, 
      \eqn{\delta}{delta}, \eqn{\rho_1}{rho_1}, \eqn{\rho_2}{rho_2}, 
      \eqn{\rho}{rho}, and eventually \eqn{\nu}{nu}
      in case of the nested CES with 4 explanatory variables,
      where in all cases the \eqn{\nu}{nu} is only required if the model
      has variable returns to scale.
      If the vector is named, the names must be \code{"gamma"},
      \code{"delta"}, \code{"rho"}, and eventually \code{"nu"}
      in case of two expanatory variables,
      \code{"gamma"}, \code{"delta_1"}, \ldots, \code{"delta_N"},
      \code{"rho"}, and eventually \code{"nu"}
      in case of the non-nested CES with \eqn{N>2} explanatory variables,
      and \code{"gamma"}, \code{"delta_1"}, \code{"delta_2"},
      \code{"rho_1"}, \code{"rho_2"}, \code{"rho"}, and eventually \code{"nu"}
      in case of the nested CES with 4 explanatory variables,
      where the order is irrelevant in all cases.}
   \item{tName}{optional character string specifying the name of the 
      time variable (\eqn{t}).}
   \item{nested}{logical. ;
      if \code{FALSE} (the default), the original CES for \eqn{n} inputs 
      proposed by Kmenta (1967) is used;
      if \code{TRUE}, the nested version of the CES
      for 3 or 4 inputs proposed by Sato (1967) is used.} 
   \item{rhoApprox}{if the absolute value of the coefficient \eqn{\rho}{rho},
      \eqn{\rho_1}{rho_1}, or \eqn{\rho_2}{rho_2} 
      is smaller than or equal to this argument,
      the endogenous variable is calculated
      using a first-order Taylor series approximation 
      at the point \eqn{\rho}{rho} = 0 (for non-nested CES functions)
      or a linear interpolation (for nested CES functions),
      because this avoids large numerical inaccuracies
      that frequently occur in calculations with non-linear CES functions
      if \eqn{\rho}{rho}, \eqn{\rho_1}{rho_1}, or \eqn{\rho_2}{rho_2}
      have very small values (in absolute terms).}
}

\value{
   A numeric vector with length equal to the number of rows of the data set
   specified in argument \code{data}.
}

\references{
Kmenta, J. (1967):
On Estimation of the CES Production Function.
International Economic Review 8, p. 180-189.

Sato, K. (1967): 
A Two-Level Constant-Elasticity-of-Substitution Production Function.
Review of Economic Studies 43, p. 201-218.
}

\seealso{\code{\link{cesEst}}.}

\author{Arne Henningsen and Geraldine Henningsen}

\examples{
   data( germanFarms, package = "micEcon" )
   # output quantity:
   germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
   # quantity of intermediate inputs
   germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput


   ## Estimate CES: Land & Labor with fixed returns to scale
   cesLandLabor <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms )

   ## Calculate fitted values
   cesCalc( c( "land", "qLabor" ), germanFarms, coef( cesLandLabor ) )


   # variable returns to scale
   cesLandLaborVrs <- cesEst( "qOutput", c( "land", "qLabor" ), germanFarms,
      vrs = TRUE )

   ## Calculate fitted values
   cesCalc( c( "land", "qLabor" ), germanFarms, coef( cesLandLaborVrs ) )
}

\keyword{models}
