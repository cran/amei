\name{optvac}
\alias{optvac}

\title{ Optimal (Static) Vaccination Policy }
\description{
  A Monte Carlo method is used to calculate the expected
  costs of a range of static vaccination policies for an
  epidemic with a known parameterization and initialization
}
\usage{
optvac(init, params, vacgrid, costs, T = 40, MCvits = 100,
     midepi = FALSE, start = 7)
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
  mortality probability, respectively }
\item{vacgrid}{ a \code{list} containing vector entries \code{$fracs} and
  \code{$stops} indicating the permissible fractions (in [0,1])
  of the population to be vaccinated and the (positive integer) of
  stopping thresholds having a maximum of \code{init$S0} }
\item{costs}{ a \code{list} containing scalar entries \code{$vac},
  \code{$death}, and \code{$infect} depicting the costs associated
  with a single vaccination, death, or the daily cost of maintaining
  an infected individual, respectively}
  \item{T}{ the maximum number of time steps during which the epidemic
  is allowed to evolve }
\item{MCvits}{ the number of Monte Carlo iterations of forward epidemic
  evolution used to determine the optimal vaccination
  policy }
\item{midepi}{ a debugging \code{logical} indicating whether to
  signal that a trajectory is unlikely to be an epidemic }
\item{start}{ at what time, after time 1 where the state is given
  by \code{init}, should vaccinations be allowed to start }
}
\details{
  This function use a Monte Carlo experiment to calculate the expected
  costs over a range of vaccination policies specified by permissible
  fractions of individuals to be vaccinated and stopping thresholds.
  These policies are constructed by simulating SIR-modeled epidemics
  that evolve according to \code{params} starting in the
  \code{init}ial configuration provided.
  The output is an object of class \code{"optvac"}, so the cost grid,
  or matrix, can be visualized with the \code{plot.optvac} generic method.
  The \code{getpolicy} function can be used to select out the best
  (and worst) one(s).

  For more details on the parameterization and simulation of the
  SIR model, etc., see \code{vignette("amei")}
}
\value{
  \code{optvac} returns an object of class \code{"optvac"}, which is a
  \code{list} containing the following components.

  \item{vacgrid }{ a copy of the input \code{vacgrid} }
  \item{C }{ a matrix of expected costs estimated for each combination
    of \code{vacgrid$fracs} and \code{vacgrid$stops} }
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
Daniel Merl <dan@stat.duke.edu>\cr
Leah R. Johnson <leah@statslab.cam.ac.uk>\cr
Robert B. Gramacy <bobby@statslab.cam.ac.uk>\cr
and Mark S. Mangel <msmangl@ams.ucsc.edu>
}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{plot.optvac}}, \code{\link{getpolicy}}, \code{\link{MCepi}} }
\examples{
## same inputs as in the MCepi example
truth <- list(b=0.00218, k=10, nu=0.4, mu=0) 
init <- list(S0=762, I0=1, R0=0, D0=0) 
costs <- list(vac=2, death=4, infect=1)

## construct a grid of valid vaccination strategies
vacgrid <- list(fracs=seq(0,1.0,0.1), stops=seq(2,init$S0-50,50))

## calculate the cost surface of all combinations in vacgrid
out.optvac <- optvac(init, truth, vacgrid, costs)

## extract the best and worst (static) policy
best <- getpolicy(out.optvac)
worst <-  getpolicy(out.optvac, which="worst")
rbind(best, worst)

## plot the cost surface along with the best and worst policy
plot(out.optvac)

## now return to MCepi for a cost comparison to no vaccination
## using these values
vac.opt <- best[3:4]
vac.opt
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ iteration }
\keyword{ optimize }
