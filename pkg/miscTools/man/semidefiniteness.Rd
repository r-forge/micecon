\name{isSemidefinite}
\alias{isSemidefinite}
\alias{isSemidefinite.default}
\alias{isSemidefinite.list}
\alias{isSemidefinite.matrix}
\alias{semidefiniteness}
\title{Positive or Negative Semidefiniteness}

\description{
   Check whether a symmetric matrix is positive or negative semidefinite.
}

\usage{
isSemidefinite( m, \dots )

\method{isSemidefinite}{default}( m, \dots )

\method{isSemidefinite}{matrix}( m, positive = TRUE,
   tol = 100 * .Machine$double.eps, method = "det", \dots )

\method{isSemidefinite}{list}( m, \dots )

semidefiniteness( m, \dots )
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
   \item{\dots}{further arguments of \code{isSemidefinite.list}
      are passed to \code{isSemidefinite.matrix};.
      further arguments of \code{semidefiniteness}
      are passed to \code{isSemidefinite};
      further arguments of other functions are currently ignored.}
}

\details{
   Function \code{semidefiniteness()} passes all its arguments
   to \code{isSemidefinite()}.
   It is only kept for backward-compatibility
   and may be removed in the future.

   If argument \code{positive} is set to \code{FALSE},
   \code{isSemidefinite()} checks for negative semidefiniteness
   by checking for positive semidefiniteness
   of the negative of argument \code{m}, i.e. \code{-m}.

   If method \code{"det"} is used (the default),
   \code{isSemidefinite()} checks whether all principal minors
   (not only the leading principal minors)
   of the matrix \code{m}
   (or of the matrix \code{-m} if argument \code{positive} is \code{FALSE})
   are larger than \code{-tol}.
   Due to rounding errors,
   which are unavoidable on digital computers,
   the calculated determinants of singular (sub-)matrices
   (which should theoretically be zero)
   can considerably deviate from zero.
   In order to reduce the probability of incorrect results
   due to rounding errors,
   \code{isSemidefinite()} does not calculate the determinants
   of (sub-)matrices with reciprocal condition numbers
   smaller than argument \code{tol}
   but sets the corresponding principal minors to (exactly) zero.

   If method \code{"eigen"} is used,
   \code{isSemidefinite()} checks whether all eigenvalues
   of the matrix \code{m}
   (or of the matrix \code{-m} if argument \code{positive} is \code{FALSE})
   are larger than \code{-tol}.
   Due to rounding errors,
   which are unavoidable on digital computers,
   those eigenvalues of a singular matrix
   that should theoretically be zero
   can considerably deviate from zero.
   In order to reduce the probability of incorrect results
   due to rounding errors,
   \code{isSemidefinite()} does not calculate the eigenvalues
   of an NxN matrix with reciprocal condition number
   smaller than argument \code{tol}
   but checks whether all N (N-1)x(N-1) submatrices
   with row i and column i, i = 1, ..., N, removed
   are positive semidefinite.
   If necessary, this procedure is done recursively.

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
