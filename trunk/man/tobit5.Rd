\name{tobit5}
\alias{tobit5}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tobit5(selection, formula2, formula3, data = sys.frame(sys.parent()), beta0 = NULL, print.level = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{selection}{ ~~Describe \code{selection} here~~ }
  \item{formula2}{ ~~Describe \code{formula2} here~~ }
  \item{formula3}{ ~~Describe \code{formula3} here~~ }
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{beta0}{ ~~Describe \code{beta0} here~~ }
  \item{print.level}{ ~~Describe \code{print.level} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
  ~~ If necessary, more details than the __description__  above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ }

 ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{
   ## uniform observable covariates (note that exclusion restriction is
   ## fulfilled)!
   ## The errors are uncorrelated, so we could use simple OLS instead.
   N <- 1000
   z <- runif(N)
   x <- runif(N)
   ## latent selection variable
   y1s <- 2*z + rnorm(N)
   y1 <- y1s > 0
   ## latent regressions
   y2s <- x + 0.5*rnorm(N)
   y3s <- -x + 0.2*rnorm(N)
   ## instead of three rnorm(), we should use rmvnorm for correlated errors.
   y <- ifelse(y1==0, y2s, y3s)
   plot(z[y1==0], y[y1==0])
   points(z[y1==1], y[y1==1], col=2)
   a <- tobit5(y1~z, y~x, y~x, print.level=print.level)
   summary(a)
}
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
