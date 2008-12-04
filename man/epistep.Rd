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
epistep(SIR, true = list(b = 0.00218, k = 10, nu = 0.4, mu = 0))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SIR}{ a \code{list} with the current scalar values of the
    number of susceptibles (\code{$S}), infecteds (\code{$I}) and
    recovereds (\code{$R}) }
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

  The epidemic described by the default parameterization ({\tt true})
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
\emph{A statistical framework for the adaptive management of 
epidemiological interventions} (2008). 
Daniel Merl, Leah R. Johnson, Robert B. Gramacy, and Marc S. Mangel.
Duke Working Paper 08-29.
\url{http://ftp.stat.duke.edu/WorkingPapers/08-29.html}

Murray, J. D. (2002) \emph{Mathematical Biology I: An Introduction}.
Springer Verlag
}
\author{
Daniel Merl <dan@stat.duke.edu>, 
Leah R. Johnson <leah@statslab.cam.ac.uk>,
Robert B. Gramacy <bobby@statslab.cam.ac.uk>,
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
