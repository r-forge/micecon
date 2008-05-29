\name{translogCheckMono}
\alias{translogCheckMono}
\title{Monotonicity of a Translog Function}

\description{
   Check monotonicity of a translog function.
}

\usage{
translogCheckMono( xNames, data, coef, increasing = TRUE,
   strict = FALSE, quadHalf = TRUE, dataLogged = FALSE )
}

\arguments{
   \item{xNames}{a vector of strings containing the names of the
      independent variables.}
   \item{data}{dataframe containing the data.}
   \item{coef}{vector containing all coefficients.}
   \item{increasing}{logical. Check whether the function is monotonically
      increasing (default, TRUE) or decreasing (FALSE)?}
   \item{strict}{logical. Check for strict (TRUE) or non-strict
      (default, FALSE) monotonicity?}
   \item{quadHalf}{logical. Multiply the quadratic terms by one half?}
   \item{dataLogged}{logical. Are the values in \code{data} already logged?}
}

\value{
   a list of class \code{translogCheckMono} containing following objects:
   \item{obs}{a vector indicating whether monotonicity is fulfilled
      at each observation.}
   \item{exog}{data frame indicating whether monotonicity is fulfilled
      for each exogenous variable at each observation.}
}

\seealso{
   \code{\link{translogEst}} and \code{\link{translogDeriv}}
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\examples{
   data( germanFarms )
   # output quantity:
   germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
   # quantity of variable inputs
   germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
   # a time trend to account for technical progress:
   germanFarms$time <- c(1:20)

   # estimate a translog production function
   estResult <- translogEst( "qOutput", c( "qLabor", "land", "qVarInput", "time" ),
      germanFarms )

   # check whether the production function is monotonically increasing
   translogCheckMono( c( "qLabor", "land", "qVarInput", "time" ),
      germanFarms, coef( estResult ) )
}

\keyword{models}
