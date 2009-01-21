\name{quadFuncEla}
\alias{quadFuncEla}
\title{Elasticities of a Quadratic Function}

\description{
   Calculate elasticities of a quadratic function.
}

\usage{
quadFuncEla( xNames, data, coef, yName = NULL, 
   shifterNames = NULL, homWeights = NULL, quadHalf = TRUE )
}

\arguments{
   \item{xNames}{a vector of strings containing the names of the
      independent variables.}
   \item{data}{dataframe containing the data.}
   \item{coef}{vector containing all coefficients.}
   \item{yName}{an optional string containing the name of the dependent
      variable.
      If it is \code{NULL}, the dependent variable is calculated
      from the independent variables and the coefficients.}
   \item{shifterNames}{an optional vector of strings containing the names
      of the independent variables that are included
      as shifters only (not in quadratic or interaction terms).}
   \item{homWeights}{named numeric vector with weighting factors
      for calculating an index that is used to normalize the variables
      for imposing homogeneity of degree zero in these variables
      (see documentation of \code{\link{quadFuncEst}}).}
   \item{quadHalf}{logical. Multiply the quadratic terms by one half?}
}

\value{
   A data.frame of class \code{quadFuncEla}, 
   where each column corresponds to one of the independent variables.
}

\seealso{
   \code{\link{quadFuncEst}}, \code{\link{quadFuncDeriv}} nad
   \code{\link{quadFuncCalc}}
}

\author{Arne Henningsen}

\examples{
   data( germanFarms )
   # output quantity:
   germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
   # quantity of variable inputs
   germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
   # a time trend to account for technical progress:
   germanFarms$time <- c(1:20)

   # estimate a translog production function
   estResult <- quadFuncEst( yName = "qOutput", 
      xNames = c( "qLabor", "land", "qVarInput", "time" ),
      data = germanFarms )

   # compute the partial production elasticities with "fitted" output
   elaFit <- quadFuncEla( xNames = c( "qLabor", "land", "qVarInput", "time" ),
      data = germanFarms, coef = coef( estResult ) )
   elaFit

   # compute the partial production elasticities with observed output
   elaObs <- quadFuncEla( xNames = c( "qLabor", "land", "qVarInput", "time" ),
      data = germanFarms, coef = coef( estResult ), yName = "qOutput" )
   elaObs
}

\keyword{models}
