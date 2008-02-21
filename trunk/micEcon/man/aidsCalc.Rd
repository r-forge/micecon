\name{aidsCalc}
\alias{aidsCalc}
\title{ Shares and Quantities of the Almost Ideal Demand System }

\description{
  Given prices, total expenditure and coefficients this function
  calculates the demanded quantities and expenditure shares based
  on the Almost Ideal Demand System.
}

\usage{
   aidsCalc( priceNames, totExpName, coef, data = NULL, priceIndex = "TL", lnp = NULL )
}

\arguments{
   \item{priceNames}{a vector of strings containing the names of the prices.}
   \item{totExpName}{a string containing the variable name of total expenditure.}
   \item{coef}{a list containing the coefficients alpha, beta, gamma,
      and (only for the translog price index) alpha0.}
   \item{data}{a data frame containing the data.}
   \item{priceIndex}{prices index (see \code{\link{aidsPx}}).}
   \item{lnp}{log of the price index (only required if '\code{priceIndex}'
      is neither 'TL' or 'S').}
}

\value{
   \code{aidsCalc} returns a list with following elements:
   \item{shares}{a data frame containing the calculated expenditure shares.}
   \item{quantities}{a data frame containing the calculated quantites.}
}

\references{
   Deaton, A.S. and J. Muellbauer (1980)
   An Almost Ideal Demand System.
   \emph{American Economic Review}, 70, p. 312-326.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{ \code{\link{aidsEst}}, \code{\link{aidsPx}} }

\examples{
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   priceNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
   shareNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )

   ## LA-AIDS
   estResult <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, priceIndex = "S" )

   # using observed shares in the Stone index
   lnp <- aidsPx( "S", priceNames, shareNames, Blanciforti86 )
   fitted <- aidsCalc( priceNames, "xFood", coef = estResult$coef,
      data = Blanciforti86, lnp = lnp )
   fitted$shares  # equal to estResult$wFitted
   fitted$quant   # equal to estResult$qFitted

   # using fitted shares in the Stone index
   fitted <- aidsCalc( priceNames, "xFood", coef = estResult$coef,
      data = Blanciforti86, priceIndex = "S" )

   ## AIDS
   estResult <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, method = "IL" )

   fitted <- aidsCalc( priceNames, "xFood", coef = estResult$coef,
      data = Blanciforti86 )

   fitted$shares  # equal to estResult$wFitted
   fitted$quant   # equal to estResult$qFitted

}

\keyword{ models }
