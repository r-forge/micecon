\name{cesCalc}
\alias{cesCalc}
\title{Calculate CES function}

\description{
   Calculate the endogenous variable
   of a \sQuote{Constant Elasticities of Scale} (CES)
   function with two explanatory variables

   \deqn{y = \gamma \: ( \delta \: x_1^{-\rho} + ( 1 - \delta ) \: x_2^{-\rho} )
      ^{-\frac{\phi}{\rho}}}{%
      y = gamma * ( delta * x1^(-rho) + ( 1 - delta ) * x2^(-rho) )^(-phi/rho)}

   and with \eqn{N} explanatory variables

   \deqn{y = \gamma \: \left( \sum_{i=1}^N \delta_i \: x_i^{-\rho} \right)
      ^{-\frac{\phi}{\rho}}}{%
      y = gamma * ( sum(i=1 to N) delta_i * x_i^(-rho) )^(-phi/rho)}

   where in the latter case
   \eqn{\sum_{i=1}^N \delta_i = 1}{sum(i=1 to N) delta_i = 1}.

   In both cases, the elesticity of substitution is
   \eqn{s = 1 / ( 1 + \rho )}{s = 1 / ( 1 + rho )}

}

\usage{
cesCalc( xNames, data, coef )
}

\arguments{
   \item{xNames}{a vector of two strings containing the names of the
      explanatory variables.}
   \item{data}{data frame containing the explanatory variables.}
   \item{coef}{numeric vector containing the coefficients of the CES:
      if the vector is unnamed,
      the order of the coefficients must be
      \eqn{\gamma}{gamma}, \eqn{\delta}{delta}, \eqn{\rho}{rho},
      and eventually \eqn{\phi}{phi}
      in case of two expanatory variables,
      \eqn{\gamma}{gamma}, \eqn{\delta_1}{delta_1}, \ldots,
      \eqn{\delta_N}{delta_N}, \eqn{\rho}{rho}, and eventually \eqn{\phi}{phi}
      in case of \eqn{N>2} explanatory variables,
      where the \eqn{\phi}{phi} is only required if the model
      has variable returns to scale;
      if the vector is named, the names must be \code{"gamma"},
      \code{"delta"}, \code{"rho"}, and eventually \code{"phi"}
      in case of two expanatory variables,
      and \code{"gamma"}, \code{"delta_1"}, \ldots, \code{"delta_N"},
      \code{"rho"}, and eventually \code{"phi"}
      in case of \eqn{N>2} explanatory variables,
      where the order is irrelevant in both cases.}
}

\value{
   A numeric vector with length equal to the number of rows of the data set
   specified in argument \code{data}.
}

\seealso{\code{\link{cesEst}}.}

\author{Arne Henningsen}

\examples{
   data( germanFarms )
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
