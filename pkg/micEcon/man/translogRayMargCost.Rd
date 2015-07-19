\name{translogRayMargCost}
\alias{translogRayMargCost}
\title{Marginal Costs of Translog Ray Production Function}

\description{
   Calculate the marginal costs of the outputs
   from a multiple-output translog ray production function.
}

\usage{
translogRayMargCost( yNames, xNames, wNames, data, coef,
   dataLogged = FALSE )
}

\arguments{
   \item{yNames}{a vector of two character strings
      containing the names of the output quantities.}
   \item{xNames}{a vector of strings containing the names of the
      input quantities.}
   \item{wNames}{a vector of strings containing the names of the
      input prices.}
   \item{data}{dataframe containing the data.}
   \item{coef}{vector containing all coefficients:
      if there are \code{n} inputs in \code{xNames},
      the \code{n+1} alpha coefficients must have names
      \code{a_0}, \ldots, \code{a_n}
      and the \code{n*(n+1)/2} beta coefficients must have names
      \code{b_1_1}, \ldots, \code{b_1_n}, \ldots, \code{b_n_n}
      (only the elements of the upper right triangle of the beta matrix
      are directly obtained from \code{coef};
      the elements of the lower left triangle are obtained by assuming
      symmetry of the beta matrix).}
   \item{dataLogged}{logical. Are the values in \code{data} already logged?}
}

\value{
   A data frame containing the marginal costs of producing the outputs.
}

\seealso{
   \code{\link{translogRayEst}}, \code{\link{translogProdFuncMargCost}},
   and \code{\link{translogCostEst}}.
}

\author{Arne Henningsen and Geraldine Henningsen}

\examples{
   data( germanFarms )
   # quantity of crop outputs:
   germanFarms$qCrop <- germanFarms$vCrop / germanFarms$pOutput
   # quantity of animal outputs:
   germanFarms$qAnimal <- germanFarms$vAnimal / germanFarms$pOutput
   # quantity of variable inputs
   germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
   # a time trend to account for technical progress:
   germanFarms$time <- c(1:20)

   # generate (artificial) prices
   germanFarms$pLand <- 200 + 15 * germanFarms$time
   germanFarms$pTime <- 1

   # estimate a two-output translog ray production function
   estResultRay <- translogRayEst( yNames = c( "qCrop", "qAnimal" ),
      xNames = c( "qLabor", "land", "qVarInput" ),
      data = germanFarms )

   # compute the marginal costs of producing the two outputs
   margCostRay <- translogRayMargCost( yNames = c( "qCrop", "qAnimal" ),
      xNames = c( "qLabor", "land", "qVarInput" ),
      wNames = c( "pLabor", "pLand", "pVarInput" ),
      data = germanFarms, coef = coef( estResultRay ) )
}

\keyword{models}
