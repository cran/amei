\name{MCepi}
\alias{MCepi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Monte Carlo Epidemic Simulation with Static Vaccination }
\description{
  This function performs a Monte Carlo simulation of epidemic
  trajectories under a static vaccination strategy.
  Statistics are collected for
  a collection of  characteristics including the state
  of the epidemic (evolution of susceptibles (S), infecteds (I),
  recovereds (R), and deads) and the cost of the vaccination
  strategy employed
}
\usage{
MCepi(init, params, vac, costs, T = 40, MCreps = 1000,
      quant = c(0.025, 0.975), midepi = FALSE, start = 7, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init}{ a \code{list} containing scalar entries \code{$S0},
    \code{$I0}, \code{$R0}, \code{$D0} depicting the initial number of
  susceptibles, infecteds, recovereds, and deads in the outbreak }
  \item{params}{ a \code{list} containing scalar entries
    \code{$b}, \code{$k}, \code{$nu}, and \code{$mu} depicting the
    \emph{true} parameters of the SIR model representing the transmission
    probability, clumpiness parameter, the recovery probability, and the
    mortality probability, respectively; optionally \code{params}
    may be an \code{\link{epistep}} function}
\item{vac}{ a \code{list} containing scalar entries \code{$frac}
  and \code{$stop} depicting the fraction to be vaccinated at each
  time step, and the vaccination (stopping) threshold, respectively }
\item{costs}{ a \code{list} containing scalar entries \code{$vac},
  \code{$death}, and \code{$infect} depicting the costs associated
  with a single vaccination, death, or daily cost of maintaining an
  infected individual, respectively}
\item{T}{ the maximum number of time steps during which the epidemic
  is allowed to evolve }
\item{MCreps}{ the number of times to repeat the Monte Carlo experiment,
  each time starting with the state in \code{init} and collecting
  characteristics of the epidemic trajectories and (vaccination/death)
  costs }
\item{quant}{ a 2-vector of quantiles to use in order to capture
  the spread in the density of characteristics of the epidemic
  trajectory and costs }
\item{midepi}{ a debugging \code{logical} indicating whether to
  signal that a trajectory is unlikely to be an epidemic }
\item{start}{ at what time, after time 1 where the state is given
  by \code{init}, should vaccinations be allowed to start }
\item{...}{ additional arguments passed to a user-defined
  \code{\link{epistep}} function that may be passed in via \code{params}}
}
\details{
  This function simulates many (\code{MCreps}) trajectories of an
  epidemic starting out in a particular state (\code{init}) and
  evolving according to a particular (\code{true}) parameterization
  under a fixed vaccination strategy (\code{vac}).  It returns
  a summary of characteristics of the state trajector(y/ies)
  and the associated \code{cost}s.  The output can be
  visualized with the generic \code{\link{plot.MCepi}} method
  and costs can be extracted with \code{\link{getcost}}.

  For more details on the parameterization and simulation of the
  SIR model, etc., and the calculation of the optimal vaccination
  strategy, please see \code{vignette("amei")}
}
\value{
  \code{MCepi} returns an object of class \code{"MCepi"}, which is a
  \code{list} containing the components listed below.

  \item{Q1 }{ a \code{data.frame} containing 6 columns (\code{S},
    \code{I}, \code{R}, \code{D}, \code{V}, \code{C}) depicting
    the first quantile (\code{quant[1]}) of the distribution of
    the evolution of the state of the epidemic (SIRD), the number
    of vaccinations (V), and the cost (C), at each time point }
  \item{Mean }{ same as \code{Q1} except the mean rather than
  a quantile }
\item{Q3}{same as \code{Q1} except the third quantile (\code{quant[2]}) }

These quantities can be visually inspected using the
\code{\link{plot.MCepi}} method 
}
\references{
\emph{A statistical framework for the adaptive management of 
epidemiological interventions} (2008). 
Daniel Merl, Leah R. Johnson, Robert B. Gramacy, and Marc S. Mangel.
Duke Working Paper 08-29.
\url{http://ftp.stat.duke.edu/WorkingPapers/08-29.html}
}
\author{
Daniel Merl <dan@stat.duke.edu>, 
Leah R. Johnson <leah@statslab.cam.ac.uk>, 
Robert B. Gramacy <bobby@statslab.cam.ac.uk>,
and Mark S. Mangel <msmangl@ams.ucsc.edu>
}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{plot.MCepi}}, \code{\link{optvac}} }
\examples{
## true epidemic parameters, initial values, and
## vaccination costs
truth <- list(b=0.00218, k=10, nu=0.4, mu=0) 
init <- list(S0=762, I0=1, R0=0, D0=0) 
costs <- list(vac=2, death=4, infect=1)

## trivial vaccination strategy -- dont vaccinate
vac <- list(frac=0, stop=0)

## simulate the resulting trajectories
init.MCepi <- MCepi(init, truth, vac, costs)

## plot the distribution of trajectories and costs
## under no vaccination
plot(init.MCepi, main="no vaccination")
plot(init.MCepi, type="costs")

## Now try the optimal strategy.
## See the optvac function for more info
vac.opt <- list(frac=0.7, stop=502)
opt.MCepi <- MCepi(init, truth, vac.opt, costs)

## plot the distribution of trajectories and costs
## under the optimal (static) vaccination
plot(opt.MCepi, main="optimal static vaccination")
plot(opt.MCepi, type="costs")

## show the total number of vaccinations
## median and quantiles
getvac(opt.MCepi)

## compare the median costs of the the initial
## (no vaccination) strategy versus the optimal
## (static) policy
T <- length(opt.MCepi$Median$C)
optC <- getcost(opt.MCepi)
initC <- getcost(init.MCepi)
rbind(initC,optC)
}
\keyword{ methods }
\keyword{ iteration }% __ONLY ONE__ keyword per line
