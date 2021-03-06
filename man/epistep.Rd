\name{epistep}
\alias{epistep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Evolve One Step of an Epidemic in Time}
\description{
  This function takes the current state of an epidemic,
  described by the values of SIR, and evolves the epidemic
  by one time step, stochastically, according to the parameterization
  provided
}
\usage{
epistep(SIR, last = NULL, true = list(b = 0.00218, k = 10, nu = 0.4, mu = 0))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SIR}{ a \code{list} with the current scalar values of the
    number of susceptibles (\code{$S}), infecteds (\code{$I}) and
    recovereds (\code{$R}) }
  \item{last}{ a (dummy) argument used to pass additional information
    necessary for the sampling the dynamics of the epidemic;
    the return-value of a user-defined \code{epistep} function would
    be automatically passed in here by the \code{\link{manage}}
    function.  This is not required by the simple SIR model
    implemented here in the default version}
  \item{true}{ a \code{list} containing scalar entries indicating
    the \emph{true} parameters according to which the SIR model evolves
    stochastically: \code{$b}, \code{$k}, \code{$nu}, and \code{$mu}
    representing the transmission probability, clumpiness parameter,
    the recovery probability, and the mortality probability, respectively }
}
\details{
  This function is intended to be passed as an argument to the
  \code{\link{manage}} function, to describe the default evolution
  of an epidemic under the SIR model.  Other, user-defined, functions
  undergoing different disease dynamics should follow the protocol (i.e.,
  inputs and outputs) prototyped by this function.  Similarly, this
  function may be used as input to \code{\link{MCmanage}} which
  depends on the \code{\link{manage}} function.

  The epidemic described by the default parameterization
  (given by \code{true})
  is an approximation of an influenza epidemic in a British
  boarding school described by Murray (see references below).
  
  For more details on the parameterization and simulation of the
  SIR model, etc., see \code{vignette("amei")}
}
\value{
  \code{epistep} returns a \code{list} containing the
  scalar integer components listed below indicating the number of
  individuals which are

  \item{rem }{ newly removed }
  \item{rec }{ newly recovered }
  \item{infect }{ newly infected }
  \item{dead }{ newly dead }
}
\references{
D. Merl, L.R. Johnson, R.B. Gramacy, and M.S. Mangel (2010).
\dQuote{\pkg{amei}: An \R Package for the Adaptive Management of
Epidemiological Interventions}.  \emph{Journal of Statistical Software}
\bold{36}(6), 1-32. \url{http://www.jstatsoft.org/v36/i06/}
  
D. Merl, L.R. Johnson, R.B. Gramacy, M.S. Mangel (2009). \dQuote{A
Statistical Framework for the Adaptive Management of Epidemiological
Interventions}.  \emph{PLoS ONE}, \bold{4}(6), e5807.
\url{http://www.plosone.org/article/info:doi/10.1371/journal.pone.0005807}

Murray, J. D. (2002) \emph{Mathematical Biology I: An Introduction}.
Springer Verlag
}
\author{
Daniel Merl <danmerl@gmail.com>\cr
Leah R. Johnson <lrjohnson@uchicago.edu>\cr
Robert B. Gramacy <rbgramacy@chicagobooth.edu>\cr
and Mark S. Mangel <msmangl@ams.ucsc.edu>
}
\seealso{ \code{\link{manage}}, \code{\link{MCmanage}} }
\examples{
## parameters to epistep (similar default except mu != 0)
true <- list(b = 0.00218, k = 0.1, nu = 0.4, mu = 0.1)
SIR <- list(S=700, I=200, R=100)

## examine the distribution of the outputs of epistep
T <- 1000
na <- rep(NA, T)
out <- data.frame(rem=na, rec=na, infect=na, dead=na)
for(t in 1:T) {
  out[t,] <- epistep(SIR=SIR, true=true)
}

## make histograms of the output
par(mfrow=c(2,2))
hist(out$rem)
hist(out$rec)
hist(out$infect)
hist(out$dead)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
