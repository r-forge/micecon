\name{outlierHadi}
\alias{outlierHadi}
\title{Multivariate Outlier Detection}

\description{
   Detection of outliers in multivariate samples using the method of Hadi (1994).
}

\usage{outlierHadi( data )}

\arguments{
   \item{data}{a data frame.}
}

\details{
   All variables in the data frame that is provided in argument \code{data}
   are used for the detection of outliers.
   Since this methodology works for numeric variables only,
   all variables in the provided data frame must be numeric.
}

\value{
   A list of class \code{outlierHadi} containing following objects:
   \item{Outliers}{A vector indicating the observations (row numbers)
      that are detected as outliers}
   \item{Xbar}{Xbar.}
   \item{Covariance}{Covariance.}
   \item{mah.out}{mah.out.}
   \item{Cb}{Cb.}
   \item{Sb}{Sb.}
   \item{Distances}{Distances.}
}

\references{
Hadi, Ali S. (1994),
\emph{A Modification of a Method for the Detection of Outliers
in Multivariate Samples},
Journal of the Royal Statistical Society (B), 2, 393-396.
}

\author{S-Plus code by Ali S. Hadi, ported to R by Arne Henningsen}

\examples{
   # load data set of average German farms
   data( "germanFarms" )

   # test based on all variables (except for the year)
   outlier <- outlierHadi( germanFarms[ , -1 ] )
   outlier$Outliers
   # subset of the data without outliers
   gfNoOutlier <- germanFarms[ - outlier$Outliers, ]

   # test based on output and input values or quantities
   outlier2 <- outlierHadi( germanFarms[ ,
      c( "vOutput", "vVarInput", "qLabor", "land" ) ] )
   outlier2$Outliers
}

\keyword{array}
