\name{MCmanage}
\alias{MCmanage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Monte Carlo Epidemic Simulation with Adaptive Vaccination }
\description{
  This function performs a Monte Carlo simulation of epidemic
  trajectories under an adaptive vaccination strategy as implemented
  by the \code{manage} function.  Statistics are tallied for
  a collection of  characteristics including the state
  of the epidemic (evolution of susceptibles (S), infecteds (I),
  recovereds (R), and deads) and the cost of the vaccination
  strategy employed
}
\usage{
MCmanage(init, epistep, vacgrid, costs, 
     pinit = list(b = 0.1, k = 0.02, nu = 0.2, mu = 0.1), 
     hyper = list(bh = c(1,3), kh = c(1,3), nuh = c(1,1), muh = c(1,1)), 
     vac0 = list(frac = 0, stop = 0), T = 40, MCreps = 30, 
     MCvits = 50, MCMCpits = 1000, bkrate = 1, vacsamps = 50,
     quant = c(0.025, 0.975), start = 7, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init}{ a \code{list} containing scalar entries \code{$S0},
    \code{$I0}, \code{$R0}, \code{$D0} depicting the initial number of
    susceptibles, infecteds, recovereds, and deads in the outbreak }
  \item{epistep}{ a function which moves the epidemic ahead one
    time-step;  see \code{\link{epistep}} }
  \item{vacgrid}{ a \code{list} containing vector entries \code{$fracs},
  \code{$stops} indicating the permissible fractions (in [0,1]) of the
population to be vaccinated and the (positive integer) of
  stopping thresholds having a maximum of \code{init$S0} }
  \item{costs}{ a \code{list} containing scalar entries \code{$vac},
  \code{$death}, and \code{$infect}, depicting the costs associated
  with a single vaccination, death, or the daily cost of maintaining
  an infected individual, respectively}
\item{pinit}{ a \code{list} containing scalar entries
    \code{$b}, \code{$k}, \code{$nu}, and \code{$mu} depicting the
    \emph{initial values} of parameters of the SIR model representing
    the transmission probability, clumpiness parameter,
    the recovery probability, and the mortality probability, respectively,
    which are subsequently sampled by MCMC from the posterior}
  \item{hyper}{ a \code{list} containing 2-vector entries describing
  parameters to the prior distribution of the parameters listed in
  the \code{pinit} argument.  The prior for \code{b} follows
  a gamma distribution with parameters $\code{$bh}$ where
  the shape is given by \code{bh[1]} and scale by \code{bh[2]}.
  The prior for \code{k} specified by parameters \code{$kh} is similar.
  The prior(s) for \code{nu} and
  \code{mu} are specified through \eqn{p_r}{p_r} and \eqn{p_d}{p_d},
  respectively, which follow Beta distributions and the default
  specification is uniform.  See \code{vignette("amei")} for more
  details }
\item{vac0}{ the initial (static) vaccination policy to be used
  before estimation of parameters begins (at \code{start}).
  This is a \code{list} with scalar entries \code{$frac}
  and \code{$stop} depicting the fraction to be vaccinated at each
  time step, and the vaccination (stopping) threshold, respectively }
  \item{T}{ the maximum number of time steps during which the epidemic
    is allowed to evolve }
\item{MCreps}{ number of times to repeat the Monte Carlo experiment,
  each time starting with the state in \code{init} and collecting
  characteristics of the epidemic trajectory and (vaccination/death)
  costs }
  \item{MCvits}{ scalar number of Monte Carlo iterations of forward epidemic
  evolution used at each time step in \code{\link{manage}} to determine
  the optimal vaccination policy }
\item{MCMCpits}{ scalar number of Markov chain Monte Carlo iterations used
  at each step to estimate the SIR model parameters in
  \code{\link{manage}}}
\item{bkrate}{ number of samples of \code{b} and \code{k}, relative
  to \code{mu} and \code{nu} before a sample of all four parameters
  is saved in \code{\link{manage}}}
\item{vacsamps}{ used to thin the MCMC samples of the parameters sampled
  from the posterior that are used to calculate optimal vaccination policies;
  this should be an integer scalar such that \code{0 < vacsamps <= MCMCpits} }
\item{quant}{ a 2-vector of quantiles to use in order to capture
  the spread in the density of characteristics of the epidemic
  trajectory and costs }
  \item{start}{ at what time, after time 1 where the state is given
    by \code{init}, should vaccinations be allowed to start }
  \item{...}{ additional arguments passed to a user-defined
  \code{\link{epistep}} function }
}
\details{
  This function simulates many (\code{MCreps}) trajectories of an
  epidemic starting out in a particular state (\code{init}) and
  evolving according the dynamics encoded in \code{\link{epistep}}
  (or some other user-defined function) under an adaptive
  vaccination strategy as implemented by \code{\link{manage}}.
  Many of the arguments to this function are simply passed
  to \code{\link{manage}}.

  It returns a summary of characteristics of the state trajector(y/ies)
  and the associated \code{cost}s.  The output can be
  visualized with the generic \code{\link{plot.MCepi}} method
  and costs can be extracted with \code{\link{getcost}}.

  For more details on the parameterization and simulation of the
  SIR model, etc., and the calculation of the optimal vaccination
  strategy, please see \code{vignette("amei")}
}
\value{
  \code{MCmanage} returns an object of class \code{"MCepi"}, which is a
  \code{list} containing the following components.

  \item{Q1 }{a \code{data.frame} containing 8 columns (\code{S},
    \code{I}, \code{R}, \code{D}, \code{V}, \code{C}, \code{frac},
    \code{stop}) depicting
    the first quantile (\code{quant[1]}) of the distribution of
    the evolution of the state of the epidemic (SIRD), the number
    of vaccinations (V), the cost (C), the fraction vaccinated (F),
    and the stopping threshold (S) at each time point }
  \item{Mean }{same as \code{Q1} except the mean rather than a quantile }
  \item{Median }{same as \code{Median} except the median rather than
    mean }
\item{Q3}{same as \code{Q1} except the third quantile (\code{quant[2]}) }

These quantities can be visually inspected using the
\code{\link{plot.MCepi}} method 
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
}
\author{
Daniel Merl <danmerl@gmail.com>\cr
Leah R. Johnson <lrjohnson@uchicago.edu>\cr
Robert B. Gramacy <rbgramacy@chicagobooth.edu>\cr
and Mark S. Mangel <msmangl@ams.ucsc.edu>
}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{manage}}, \code{\link{MCepi}},
  \code{\link{plot.MCepi}} }
\examples{
## adaptively manage the epidemic with the following
## initial population
init <- list(S0=762, I0=1, R0=0, D0=0) 

## construct a grid of valid vaccination strategies
## and specify costs
## using a smaller grid for faster check times; try the commented out
## setting for higher fidelity
vacgrid <- list(fracs=seq(0,1.0,0.25), stops=seq(2,init$S0-75,150))
## vacgrid <- list(fracs=seq(0,1.0,0.1), stops=seq(2,init$S0-75,75))
costs <- list(vac=2, death=4, infect=1)

## run the Monte Carlo management experiment, reducing MCMCpits
## for faster check times; try default (commented out) version
out.MCmanage <- MCmanage(init, epistep, vacgrid, costs, 
                         MCMCpits=100, MCreps=5)
## out.MCmanage <- MCmanage(init, epistep, vacgrid, costs)

## plot the trajectories of SIR and the associated costs
plot(out.MCmanage, main="optimal adaptive vaccination")
plot(out.MCmanage, type="costs")

## extract the distribution of the number of
## cumulative vaccinations via median and quantiles
getvac(out.MCmanage)

## plot the distribution fractions vaccinated and
## stopping times
plot(out.MCmanage, type="fracs")
plot(out.MCmanage, type="stops")

## get the final median cost and quantiles --
## these can be compared with the static ones
## calculated by MCepi
getcost(out.MCmanage)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ iteration }
\keyword{ optimize }
