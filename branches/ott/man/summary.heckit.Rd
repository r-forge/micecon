\name{summary.heckit}
\alias{summary.heckit}
\alias{summary.heckit5}
\alias{print.summary.heckit}
\alias{print.summary.heckit5}

\title{Print output of heckit estimation}

\description{
   This function returns a summary of the 2-step-Heckmann/heckit or
   tobit-5 (switching regression) estimation.
}

\usage{
   \method{summary}{heckit}(object, ...)
   \method{summary}{heckit5}(object, ...)
   \method{print}{summary.heckit}(x,
digits=max(3, getOption("digits") - 3), signif.stars=getOption("show.signif.stars"),
                                  part="full", ...)
   \method{print}{summary.heckit5}(x,
digits=max(3, getOption("digits") - 3), signif.stars=getOption("show.signif.stars"),
                                  part="full", ...)
}

\arguments{
   \item{object, x}{an object of suitable class to be summarised.}
  \item{part}{character string: which parts of the summary to print:
    "full" for all
    the estimated parameters: probit selection, outcome
    estimates, correlation and residual variance; or "outcome" for the outcome
    results only.}
  \item{digits}{numeric, (suggested) number of significant digits}
  \item{signif.stars}{logical, whether to print the significance starts
    and corresponding legend}
   \item{...}{further arguments passed to \code{\link{summary.glm}}
      and \code{\link{summary.lm}}.}
}

\details{The inverse Mill's ratios are not reported.  Use
  \code{summary(object$twoStep[1|2])} to get the full output of the
  second step OLS estimation.

  The variance-covariance matrix of the two-step estimator is currently
  only partially implemented.
}

\value{
  Summary methods return an object of class \code{summary.heckit} or
  \code{summary.heckit5}.  Print methods return the argument invisibly. 
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
print(summary(wooldridge), part="outcome")
                             # Only the outcome equation
}

\keyword{models}



