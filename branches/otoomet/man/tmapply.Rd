\name{tmapply}
\alias{tmapply}
\title{tapply a function over matrix columns}
\description{
  tapply a function over each column of a matrix, that is to each group
  inside of columns, given by a unique combination of the levels of
  certain factors.
}
\usage{
tmapply(X, INDEX, FUN = sum, ...)
}
\arguments{
  \item{X}{matrix}
  \item{INDEX}{vector of factor, the length must equal the number of
    rows of \code{X}}
  \item{FUN}{the function to be applied}
  \item{\dots}{further argument to \code{FUN}}
}
\value{
  Matrix, number of rows equal to the number of levels in \code{INDEX},
  the number of columns equal to that of \code{X}.  The elements are
  \code{FUN} applied to elements of the current column of \code{X}
  corresponding to current level of \code{INDEX}.}
}
\author{Ott Toomet \email{siim (a) obs ee}}
\seealso{\code{\link{tapply}}, \code{\link{mapply}}}
\examples{
a <- matrix(c(1,2,3,4,10,20), 3, 2, byrow=TRUE)
b <- factor(c("a", "a", "b"))
tmapply(a, b)
tmapply(a, b, mean)
}
\keyword{manip}
\keyword{utilities}
