\name{aidsConsist}
\alias{aidsConsist}
\alias{testConsist.aidsEst}
\title{Consistency Test of the AIDS}
\description{
   Test whether the specification of the AIDS is consistent
   with microeconomic demand theory (i.e. utility maximisation).
}

\usage{
aidsConsist( priceNames, totExpName, data, coef,
   priceIndex = "TL", basePrices = NULL, baseShares = NULL,
   shareNames = NULL )

\method{testConsist}{aidsEst}( object, ... )
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
   \item{object}{an object of class \code{aidsEst}.}
   \item{...}{currently not used.}
}

\details{
   \code{testConsist.aidsEst} is a wrapper function to \code{aidsConsist}
   that extracts all relevant arguments for \code{aidsConsist}
   from \code{object}.

   If argument \code{shareNames} is specified,
   observed shares are used for the calculation of the 'C' matrices
   to check for concavity;
   if argument \code{shareNames} is \code{NULL} (e.g., not specified),
   fitted shares are used for the calculation of the 'C' matrices.
}

\value{
   a list containing following elements:
   \item{monotony}{a logical vector indicating whether the monotony
      condition is fulfilled at each observation.}
   \item{mPercent}{percent of observations where the monotony condition is
      fulfilled.}
   \item{concavity}{a logical vector indicating whether the concavity
      condition is fulfilled at each observation.}
   \item{cPercent}{percent of observations where the concavity condition is
      fulfilled.}
   \item{cMatrices}{a list of the 'C' matrices for each observation to check
      for concavity (see Deaton and Muellbauer, 1980b, p.76 ).}
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

\seealso{\code{\link{aidsEst}}, \code{\link{aidsElas}}}

\examples{
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   priceNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
   shareNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )

   estResult <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, method = "IL" )
   aidsConsist( priceNames, "xFood", data = Blanciforti86,
      coef = estResult$coef, shareNames = shareNames )

   # the same can be obtained in an easier way
   testConsist.aidsEst( estResult )
}

\keyword{models}
