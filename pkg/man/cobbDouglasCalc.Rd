\name{cobbDouglasCalc}
\alias{cobbDouglasCalc}
\title{Calculate dependent variable of a Cobb-Douglas function}

\description{
   Calculate the dependent variable of a Cobb-Douglas function.
}

\usage{
cobbDouglasCalc( xNames, data, coef, dataLogged = FALSE  )
}

\arguments{
   \item{xNames}{a vector of strings containing the names of the
      independent variables.}
   \item{data}{data frame containing the data.}
   \item{coef}{vector of named elements containing the coefficients:
      if there are \code{n} exogenous variables in \code{xNames},
      the \code{n+1} coefficients must have names
      \code{a_0}, \ldots, \code{a_n}.}
   \item{dataLogged}{logical. Are the values in \code{data} already logged?}
}

\value{
   A vector containing the endogenous variable.
   If the inputs are provided as logarithmic values
   (argument \code{dataLogged} is \code{TRUE}),
   the endogenous variable is returned as logarithm;
   non-logarithmic values are returned otherwise.
}

\seealso{\code{\link{translogCalc}}.}

\author{Arne Henningsen}

\examples{
   data( germanFarms )
   # output quantity:
   germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
   # quantity of variable inputs
   germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
   # a time trend to account for technical progress:
   germanFarms$time <- c(1:20)

   # estimate a Cobb-Douglas production function
   estResult <- translogEst( "qOutput", c( "qLabor", "land", "qVarInput", "time" ),
      germanFarms, linear = TRUE )

   cobbDouglasCalc( c( "qLabor", "land", "qVarInput", "time" ), germanFarms,
      coef( estResult )[ 1:5 ] )
   #equal to estResult$fitted
}

\keyword{models}
