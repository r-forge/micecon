\name{semidefiniteness}
\alias{semidefiniteness}
\title{Positive or Negative Semidefiniteness}

\description{
   Check whether a symmetric matrix is positive or negative semidefinite.
}

\usage{semidefiniteness( m, positive = TRUE,
   tol = 100 * .Machine$double.eps, method = "det" )}

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
   Please note that a matrix can be
   neither positive semidefinite nor negative semidefinite.
}

\value{
   \code{semidefiniteness} returns a locigal value
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
   semidefiniteness( matrix( 1, 3, 3 ))

   # a negative semidefinite matrix
   semidefiniteness( matrix(-1, 3, 3 ), positive = FALSE )

   # a matrix that is positive and negative semidefinite
   semidefiniteness( matrix( 0, 3, 3 ))
   semidefiniteness( matrix( 0, 3, 3 ), positive = FALSE )

   # a matrix that is neither positive nor negative semidefinite
   semidefiniteness( symMatrix( 1:6 ) )
   semidefiniteness( symMatrix( 1:6 ), positive = FALSE )
}

\keyword{array}
