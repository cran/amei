\name{manage}
\alias{manage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Adaptive Management of Epidemiological Interventions }
\description{
  This function adaptively manages an epidemic by learning about
  (parameters determining) its behavior online and adjusting an
  optimal vaccination strategy accordingly
}
\usage{
manage(init, epistep, vacgrid, costs, T = 40, Tstop = T, 
     pinit = list(b = 0.1, k = 0.02, nu = 0.2, mu = 0.1), 
     hyper = list(bh = c(1,3), kh = c(1,3), nuh = c(1,1), muh = c(1,1)),
     vac0=list(frac=0, stop=0), MCvits = 10, MCMCpits = 1000, 
     bkrate = 1, vacsamps = 100, start = 8, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init}{ a \code{list} containing scalar entries \code{$S0},
    \code{$I0}, \code{$R0}, \code{$D0} depicting the initial number of
    susceptibles, infecteds, recovereds, and dead individuals in the
    outbreak; an \code{"epiman"} object may also be provided in order
    to continue managing an ongoing epidemic (see Details) }
  \item{epistep}{ a function which moves the epidemic ahead one
    time-step; see \code{\link{epistep}} }
  \item{T}{ the long-term horizon used to calculate optimal
    vaccination strategies }
  \item{Tstop}{ the maximum number of time steps during which the epidemic
   is allowed to evolve; default is identical to \code{T} }
 \item{vacgrid}{ a \code{list} containing vector entries \code{$fracs} and
  \code{$stops} indicating the permissible fractions (in [0,1])
  of the population to be vaccinated and the (positive integer) of
  stopping thresholds having a maximum of \code{init$S0}; can be
  \code{NULL}, see details below }
\item{costs}{ a \code{list} containing scalar entries \code{$vac},
  \code{$death} and \code{$infect} depicting the costs associated
  with a single vaccination or death, or the daily cost of maintaining
  an infected individual, respectively; can be \code{NULL}, but only
  when \code{vacgrid = NULL} also }
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
  time step, and the vaccination (stopping) threshold, respectively.
  The default corresponds to no initial vaccination }
\item{MCvits}{ scalar number of Monte Carlo iterations of forward epidemic
  evolution used at each time step to determine the optimal
  vaccination policy }
\item{MCMCpits}{ scalar number of Markov chain Monte Carlo iterations used
  at each step to estimate the SIR model parameters }
\item{bkrate}{ number of samples of \code{b} and \code{k}, relative
  to \code{mu} and \code{nu} before a sample of all four parameters
  is saved; this acknowledges that \code{b} and \code{k} are correlated
  and thus mix slower than \code{mu} and \code{nu} }
\item{vacsamps}{ used to thin the MCMC samples of the parameters sampled
  from the posterior that are used to calculate optimal vaccination policies;
  this should be an integer scalar such that \code{0 < vacsamps <= MCMCpits} }
\item{start}{ at what time, after time 1 where the state is given
  by \code{init}, should vaccinations be allowed to start }
\item{...}{ additional arguments passed to a user-defined
  \code{\link{epistep}} function }
}
\details{
  At each time step of the epidemic -- evolving stochastically according
  to the initialization in \code{init} and progression described by
  \code{epistep} -- the parameters are inferred by sampling from their
  posterior distribution conditioned on the available data via MCMC.
  These samples are fed into the Monte Carlo method for determining
  the optimal vaccination strategy for (the remainder) of the epidemic.
  That policy is then enacted, and then time is incremented.

  Parameter estimation (alone) can be performed by specifying the
  \dQuote{null} vaccination grid \code{vacgrid = NULL},
  i.e., with \code{vac0}.

  It is possible to \dQuote{continue} managing an epidemic by specifying
  the output of \code{manage}, an \code{"epiman"}-class object, as the
  \code{init} argument and specifying a finishing time (\code{T}) which
  is greater than the finishing time used to create the object (in the
  previous call to \code{manage}); see example below.
  
  For more details on the parameterization and simulation of the
  SIR model, etc., and the calculation of the optimal vaccination
  strategy, see \code{vignette("amei")}
}
\value{
  \code{management} returns an object of class \code{"epiman"},
  which is a \code{list} containing the components listed below.
  
  \item{soln }{ a \code{data.frame} describing the evolution of the
    epidemic giving the following for each time step (in columns):
    a time index (\code{$TIME}); the total number of susceptibles (\code{\$S}),
    infecteds (\code{$I}), recovereds (\code{$R}), and deads (\code{$D});
    the number of people who became infected (\code{$itilde}), recovered
    (\code{$rtilde}), or died (\code{$dtilde}), in that particular
    time step; the total number vaccinated (\code{$V}), or culled
    (\code{$QC}, not supported in this version); and the corresponding
    cumulative cost of the epidemic (\code{$C} so far) }
  \item{prime}{ a \code{data.frame} describing the changes in susceptibles
    and infecteds due to a vaccination strategy }
  \item{vacgrid }{ a copy of the input \code{vacgrid} }
  \item{pols }{ a \code{data.frame} describing the policy enacted at
    each time step with columns giving the fraction of the
    population vaccinated (\code{$frac}), and the (stopping) threshold
    (\code{$stop}) }
  \item{vactimes }{ a scalar integer vector indicating the times at
    which vaccinations actually occurred; this coincides with changes to
    \code{soln$V} }
  \item{samp }{\code{data.frame} containing samples from the
    posterior distribution of the parameters \code{$b}, \code{$k},
    \code{$nu}, and \code{$mu} collected during the final time step }
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
\note{ 
  It may be important to plot the epidemic trajectory, with
  the generic method \code{plot.epiman}, or inspect
  the output \code{$soln[T,]}, to check that the full dynamics
  of the epidemic have played out in the number of time steps,
  \code{T}, alotted
}
\seealso{ \code{\link{epistep}}, \code{\link{MCmanage}},
  \code{\link{plot.epiman}} }
\examples{
## manage an epidemic evolving according to epistep with
## with the following initial population
init <- list(S0=762, I0=1, R0=0, D0=0) 

## construct a grid of valid vaccination strategies
## and specify costs
vacgrid <- list(fracs=seq(0,1.0,0.1), stops=seq(2,init$S0-75,75))
costs <- list(vac=2, death=4, infect=1)

## adaptively manage the epidemic
out.man <- manage(init, epistep, vacgrid, costs)

## plot the trajectories of SIR and the associated costs
plot(out.man)
plot(out.man, type="cost")
getcost(out.man)

## plot the samples from the posterior distribution of
## the parameters obtained during the last time step
true <- as.list(formals(epistep)$true)
plot(out.man, type="params", true=true)

## Bobby isnt really sure what this is plotting
plot(out.man, type="fracs")
plot(out.man, type="stops")

## managing an epidemic in two stages
out.man <- manage(init, epistep, vacgrid, costs, Tstop=10)
plot(out.man)
out.man <- manage(out.man, epistep, vacgrid, costs, Tstop=20)
plot(out.man)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ optimize }
