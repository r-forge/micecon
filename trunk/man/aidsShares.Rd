\name{aidsShares}
\alias{aidsShares}
\title{ Expenditure Shares of the Almost Ideal Demand System }

\description{
  Given prices, total expenditure and coefficients this function
  calculates the demanded expenditure shares and quantities based
  on the Almost Ideal Demand System.
}

\usage{
   aidsShares( pNames, xtName, data = NULL,
   alpha0 = 0, alpha, beta, gamma, px = "TL", lnp = NULL )
}

\arguments{
   \item{pNames}{vector containing the names of the prices. }
   \item{xtName}{name of the total expenditure variable. }
   \item{data}{data frame containing the variables. }
   \item{alpha0}{coefficient \eqn{\alpha_0}. }
   \item{alpha}{vector of \eqn{\alpha} coefficients. }
   \item{beta}{vector of \eqn{\beta} coefficients. }
   \item{gamma}{matrix of \eqn{\gamma} coefficients. }
   \item{px}{prices index (see \code{\link{aidsPx}}) }
   \item{lnp}{log of the price index.}
}

\value{
   \code{aidsShares} returns a list with following elements:
   \item{shares}{a data frame containing the calculated expenditure shares}
   \item{quantities}{a data frame containing the calculated quantites}
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
   pNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
   wNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )

   ## LA-AIDS
   result <- aidsEst( pNames, wNames, "xFood",
      data = Blanciforti86, method = "LA:L" )

   lnp <- aidsPx( "L", pNames, wNames, Blanciforti86 )

   fitted <- aidsShares( pNames, "xFood", Blanciforti86,
      alpha = result$coef$alpha, beta = result$coef$beta,
      gamma = result$coef$gamma, lnp = lnp )

   fitted$shares  # equal to result$wFitted
   fitted$quant   # equal to result$qFitted

   ## AIDS
   result <- aidsEst( pNames, wNames, "xFood",
      data = Blanciforti86, method = "MK:L" )

   fitted <- aidsShares( pNames, "xFood", Blanciforti86,
      alpha = result$coef$alpha, beta = result$coef$beta,
      gamma = result$coef$gamma )

   fitted$shares  # equal to result$wFitted
   fitted$quant   # equal to result$qFitted

}

\keyword{ models }
