\name{aidsUtility}
\alias{aidsUtility}
\title{Indirect Utility Function of the Almost Ideal Demand System}

\description{
   This function calculates the utility level
   given prices and total expenditure using the indirect utility function
   of the Almost Ideal Demand System.
}

\usage{
aidsUtility( priceNames, totExpName, coef, data )
}

\arguments{
   \item{priceNames}{a vector of strings containing the names of the prices.}
   \item{totExpName}{a string containing the variable name of total expenditure.}
   \item{coef}{a list containing the coefficients in elements
      alpha0 (scalar), alpha (vector), beta (vector), gamma (matrix),
      and possibly beta0 (scalar, if not given, it is assumed to be 1).}
   \item{data}{a data frame containing the data.}
}

\value{
   \code{aidsUtility} returns a numeric vector
   that contains the utility levels.
}

\references{
   Deaton, A.S. and J. Muellbauer (1980)
   An Almost Ideal Demand System.
   \emph{American Economic Review}, 70, p. 312-326.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{aidsEst}}, \code{\link{aidsCalc}}}

\examples{
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   priceNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
   shareNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )

   ## estimate the (non-linear) AIDS
   estResult <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, method = "IL" )

   # calculate the utility levels of each year
   utility <- aidsUtility( priceNames, "xFood", coef = coef( estResult ),
      data = Blanciforti86 )
}

\keyword{ models }
