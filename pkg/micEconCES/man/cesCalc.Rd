\name{cesCalc}
\alias{cesCalc}
\title{Calculate CES function}

\description{
   Calculate endogenous variable of a Constant-Elasticities-of-Scale (CES)
   function with two exogenous variables.

   \deqn{y = \gamma * ( \alpha * x1^\rho + ( 1 - \alpha ) * x2^\rho )^{\phi/\rho}}{%
      y = gamma * ( alpha * x1^rho + ( 1 - alpha ) * x2^rho )^(phi/rho)}
   where the elesticity of substitution is
   \eqn{s = 1 / ( 1 - \rho )}{s = 1 / ( 1 - rho )}

}

\usage{
cesCalc( xNames, data, coef )
}

\arguments{
   \item{xNames}{a vector of two strings containing the names of the
      explanatory variables.}
   \item{data}{data frame containing the explanatory variables.}
   \item{coef}{numeric vector of length 3 for constant returns to scale
      (without \eqn{\phi}{phi}, i.e. \eqn{\phi}{phi} = 1)
      or of length 4 (without \eqn{\phi}{phi});
      if the vector is unnamed, the order is assumed to be
      \eqn{\gamma}{gamma}, \eqn{\alpha}{alpha}, \eqn{\rho}{rho},
      and eventually \eqn{\phi}{phi};
      if the vector is named, the names must be \code{"gamma"},
      \code{"alpha"}, \code{"rho"}, and eventually \code{"phi"},
      where the order is irrelevant.}
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
