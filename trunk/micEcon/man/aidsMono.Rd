\name{aidsMono}
\alias{aidsMono}
\title{Monotonicity of the AIDS}
\description{
   \code{aidsMono} checks whether the expenditure function
   of an estimated Almost Ideal Demand System (AIDS)
   is monotonic increasing in prices,
   which implies that all demanded quantities and expenditure shares
   are non-negative.
}

\usage{
aidsMono( priceNames, totExpName, data, coef,
   priceIndex = "TL", basePrices = NULL, baseShares = NULL,
   shareNames = NULL )
}

\arguments{
   \item{priceNames}{a vector of strings containing the names of the prices.}
   \item{totExpName}{a string containing the variable name of total expenditure.}
   \item{data}{a data frame containing the data.}
   \item{coef}{a list containing the coefficients alpha, beta, gamma,
      and (only for the translog price index) alpha0.}
   \item{priceIndex}{a character string specifying the price index
      (see \code{\link{aidsPx}})
      or a numeric vector providing the log values of the price index.}
   \item{basePrices}{a vector specifying the base prices for the
      Paasche, Laspeyres, and Tornqvist price index.}
   \item{baseShares}{a vector specifying the base expenditure shares for the
      Laspeyres, simplified Laspeyres, and Tornqvist index.}
   \item{shareNames}{a vector of strings containing the names of the expenditure
      shares.}
}

\details{
   Internally, \code{aidsMono} passes its arguments to \code{aidsCalc}
   and then checks for each observation, wether all expenditure shares
   are non-negative.
   If argument \code{shareNames} is specified,
   the price index of the LA-AIDS is calculated using observed shares;
   if argument \code{shareNames} is \code{NULL} (e.g., not specified),
   fitted shares are used for in the price index.
}

\value{
   a list containing following elements:
   \item{monotony}{a logical vector indicating whether the monotony
      condition is fulfilled at each observation.}
   \item{mPercent}{percent of observations where the monotony condition is
      fulfilled.}
}
\references{
   Deaton, A.S. and J. Muellbauer (1980a)
   An Almost Ideal Demand System.
   \emph{American Economic Review}, 70, p. 312-326.

   Deaton, A.S. and J. Muellbauer (1980b)
   \emph{Economics and Consumer Behavior},
   Cambridge University Press, Cambridge.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{aidsEst}}, \code{\link{aidsCalc}}}

\examples{
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   priceNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
   shareNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )

   estResult <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, method = "IL" )
   aidsMono( priceNames, "xFood", data = Blanciforti86,
      coef = estResult$coef )
}

\keyword{models}
