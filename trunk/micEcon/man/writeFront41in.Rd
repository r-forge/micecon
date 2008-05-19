\name{writeFront41in}
\alias{writeFront41in}
\title{Write instruction file for Frontier 4.1}

\description{
   Write an instruction file and a data file for Tim Coelli's program
   Frontier 4.1 that performs stochastic frontier analysis.
}

\usage{
   writeFront41in( data, crossSectionName, timePeriodName,
   yName, xNames = NULL, zNames = NULL,
   translog = FALSE, quadHalf = TRUE, logData = FALSE,
   functionType = 1, modelType = 1, logDepVar = TRUE, mu = FALSE, eta = FALSE,
   insFile = "front41.ins", dtaFile = sub( "\\.ins$", ".dta", insFile ),
   outFile = sub( "\\.ins$", ".out", insFile ), startUpFile = "front41.000",
   iprint = 5, indic = 1, tol = 0.00001, tol2 = 0.001, bignum = 1.0E+16,
   step1 = 0.00001, igrid2 = 1, gridno = 0.1, maxit = 100, ite = 1 )
}

\arguments{
   \item{data}{data frame that contains the data.}
   \item{crossSectionName}{string: name of the cross section identifier.}
   \item{timePeriodName}{string: name of the time period identifier.}
   \item{yName}{string: name of the endogenous variable.}
   \item{xNames}{a vector of strings containing the names of the X variables
      (exogenous variables of the production or cost function).}
   \item{zNames}{a vector of strings containing the names of the Z variables
      (variables explaining the efficiency level).}
   \item{translog}{logical (should a translog function be estimated?) or
      a vector of strings containing the names of the variables
      that should be used to construct the quadratic and interaction terms
      of the translog function.
      If this argument is \code{TRUE}, quadratic and interaction terms of all
      X variables are added, i.e. it is set equal to argument \code{xNames}.}
   \item{quadHalf}{logical. Multiply the quadratic terms of the translog
      function by one half?}
   \item{logData}{logical. Take the logarithm of the output and input
      quantities before writing them to the data file?}
   \item{functionType}{function type: either 1 for 'production function' or
      2 for 'cost function'.}
   \item{modelType}{model type: either 1 for 'Error Components Frontier' or
      2 for 'Tech. Eff. Effects Frontier'.}
   \item{logDepVar}{logical. Is the dependent variable logged.}
   \item{mu}{logical. Should a 'mu' (if modelType = 1) or a delta0
      (if modelType = 2) be included in the estimation.}
   \item{eta}{logical. Should an 'eta' be included in the estimation
      (only if modelType = 1).}
   \item{insFile}{name of the instruction file.}
   \item{dtaFile}{name of the data file.}
   \item{outFile}{name of the output file.}
   \item{startUpFile}{name of the start-up file. If this argument is NULL,
      no start-up file is written.}
   \item{iprint}{numeric. Print info every \code{iprint} iterations;
      if this argument is 0, do not print.}
   \item{indic}{numeric. Use in unidimensional search procedure:
      indic = 2 says do not scale step length in unidimensional search;
      indic = 1 says scale (to length of last step) only if last step was smaller;
      indic = any other number says scale (to length of last step).}
   \item{tol}{numeric. Convergence tolerance (proportiona).}
   \item{tol2}{tolerance used in uni-dimensional search procedure.}
   \item{bignum}{numeric. Used to set bounds on den \& dist.}
   \item{step1}{numeric. Size of 1st step in search procedure.}
   \item{igrid2}{numeric. 1 = double accuracy, 0 = single accuracy.}
   \item{gridno}{numeric. Steps taken in single accuracy grid search on gamma.}
   \item{maxit}{numeric. Maximum number of iterations permitted}
   \item{ite}{numeric. 1 = print all efficiency estimates; 
      0 = print only the mean efficiency.}
}

\value{
   writeFront41in writes an instruction file and an data file to disk.
}

\references{
   Battese, G.E. and T. Coelli (1992), Frontier production functions,
      technical efficiency and panel data: with application to paddy
      farmers in India. \emph{Journal of Productivity Analysis}, 3, 153-169.

   Battese, G.E. and T. Coelli (1995), A model for technical inefficiency effects
      in a stochastic frontier production function for panel data.
      \emph{Empirical Economics}, 20, 325-332.

   Coelli, T. (1996) A Guide to FRONTIER Version 4.1: A Computer
      Program for Stochastic Frontier Production and Cost Function
      Estimation, CEPA Working Paper 96/08,
      \url{http://www.uq.edu.au/economics/cepa/frontier.htm},
      University of New England.
}

\seealso{\code{\link{readFront41out}}}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\examples{
   data( Coelli )
   Coelli$logOutput  <- log( Coelli$output )
   Coelli$logCapital <- log( Coelli$capital )
   Coelli$logLabour  <- log( Coelli$labour )

   writeFront41in( Coelli, "firm", "time", "logOutput",
      c( "logCapital", "logLabour" ) )

   \dontrun{
   system( "wine front41.exe" )
   }
}

\keyword{models}

