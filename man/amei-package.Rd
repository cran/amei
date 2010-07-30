\name{amei-package}
\alias{amei-package}
\alias{amei}
\docType{package}
\title{Adaptive Management of Epidemiological Interventions}
\description{
  This package provides a flexible statistical framework for generating
  optimal epidemiological interventions.  The aim is to minimize the
  total expected cost of an emerging epidemic, while simultaneously
  propagating uncertainty regarding underlying disease parameters
  through to the decision process via Bayesian posterior inference.
  The strategies produced through this framework are adaptive:
  vaccination schedules are iteratively adjusted to reflect the
  anticipated trajectory of the epidemic given the current population
  state and updated parameter estimates.
}
\author{
Daniel Merl <dan@stat.duke.edu>\cr
Leah R. Johnson <leah@statslab.cam.ac.uk>\cr
Robert B. Gramacy <bobby@statslab.cam.ac.uk>\cr
and Mark S. Mangel <msmangl@ams.ucsc.edu>

Maintainer: 
Robert B. Gramacy <bobby@statslab.cam.ac.uk>
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
\keyword{ package }
