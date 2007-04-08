\name{summary.heckit}
\alias{summary.heckit}
\title{Print output of heckit estimation}

\description{
   This function returns a summary of the 2-step-Heckmann / heckit estimation.
}

\usage{
   \method{summary}{heckit}( object, part="outcome", ... )
}

\arguments{
   \item{object}{an object of class \code{heckit}.}
  \item{part}{character string: which parts to extract: "full" for all
    the estimated parameters: probit selection, outcome
    estimates and inverse Mills ratio; or "outcome" for the outcome
    results only (including the inverse Mills ratio)}
   \item{...}{further arguments passed to \code{\link{summary.glm}}
      and \code{\link{summary.lm}}.}
}

\value{
  \code{summary.heckit} returns an object of class \code{heckit}.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}, Ott
  Toomet \email{otoomet@ut.ee}}

\seealso{\code{\link{heckit}}}

\examples{
data( Mroz87 )
wooldridge <- heckit( lfp ~ nwifeinc + educ + exper + I( exper^2 ) + age + 
   kids5 + kids618, log( wage ) ~ educ + exper + I( exper^2 ), Mroz87 )
summary( wooldridge )        # summary of the 1st step probit estimation
                             # (Example 17.1, p. 562f) and 2nd step OLS
                             # regression
summary(wooldridge, part="full")
                             # full summary
}

\keyword{models}



