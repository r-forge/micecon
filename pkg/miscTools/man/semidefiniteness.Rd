\name{isSemidefinite}
\alias{isSemidefinite}
\alias{semidefiniteness}
\title{Positive or Negative Semidefiniteness}

\description{
   Check whether a symmetric matrix is positive or negative semidefinite.
}

\usage{
isSemidefinite( m, positive = TRUE,
   tol = 100 * .Machine$double.eps, method = "det" )

semidefiniteness( m, positive = TRUE,
   tol = 100 * .Machine$double.eps, method = "det" )
}

\arguments{
   \item{m}{a symmetric quadratic matrix or a list containing
      symmetric quadratic matrices.}
   \item{positive}{logical. Check for positive semidefiniteness
      (if \code{TRUE}, default)
      or for negative semidefiniteness (if \code{FALSE}).}
   \item{tol}{tolerance level (values between \code{-tol} and \code{tol}
      are considered to be zero).}
   \item{method}{method to test for semidefiniteness, either
      checking the signs of the principal minors (if \code{"det"})
      or checking the signs of the eigenvalues (if \code{"eigen"}).}
}

\details{
   Functions \code{isSemidefinite()} and \code{semidefiniteness()}
   are identical
   Function \code{semidefiniteness()} may be depreciated in the future,
   because the name of function \code{isSemidefinite()}
   seems to be more consistent to the names of other functions in R.

   If argument \code{positive} is set to \code{FALSE},
   negative semidefiniteness is checked
   by checking for potitive semidefiniteness
   of the negative of argument \code{m}, i.e. \code{-m}.

   If method \code{"det"} is used (the default),
   it is checked whether all principal minors
   (not only the leading principal minors)
   are larger than \code{-tol}.
   If the reciprocal condition number of a sub-matrix
   is smaller than argument \code{tol},
   \code{semidefiniteness()} does not calculate its determinant
   but sets the corresponding principal minor to zero.

   If method \code{"eigen"} is used,
   it is checked whether all eigenvalues
   are larger than \code{-tol}.

   Please note that a matrix can be
   neither positive semidefinite nor negative semidefinite.
}

\value{
   \code{isSemidefinite()} and \code{semidefiniteness()}
   return a locigal value
   or a logical vector (if argument \code{m} is a list)
   indicating whether the matrix (or each of the matrices)
   is positive/negative (depending on argument \code{positive})
   semidefinite.
}

\references{
   Chiang, A.C. (1984):
   \emph{Fundamental Methods of Mathematical Economics},
   3rd ed., McGraw-Hill.
   
   Gantmacher, F.R. (1959):
   \emph{The Theory of Matrices},
   Chelsea Publishing.
}

\author{Arne Henningsen}

\examples{
   # a positive semidefinite matrix
   isSemidefinite( matrix( 1, 3, 3 ))

   # a negative semidefinite matrix
   isSemidefinite( matrix(-1, 3, 3 ), positive = FALSE )

   # a matrix that is positive and negative semidefinite
   isSemidefinite( matrix( 0, 3, 3 ))
   isSemidefinite( matrix( 0, 3, 3 ), positive = FALSE )

   # a matrix that is neither positive nor negative semidefinite
   isSemidefinite( symMatrix( 1:6 ) )
   isSemidefinite( symMatrix( 1:6 ), positive = FALSE )
}

\keyword{array}
