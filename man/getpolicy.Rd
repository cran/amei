\name{getpolicy}
\alias{getpolicy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract The Optimal Static Vaccination Policy }
\description{
  Extracts the optimal static vaccination policy from an
  \code{"optvac"}-class object
}
\usage{
getpolicy(obj, which = c("best", "worst"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ and \code{"optvac"}-class object }
  \item{which}{ optionally, by supplying \code{which = "worst"}
    the worst static vaccination policy can be extracted }
}
\details{
  This function is designed to work with the output of
  the \code{\link{optvac}} function, and to provide inputs
  to the \code{\link{MCepi}} function.  It searches the cost
  grid defined by \code{vacgrid} for the lowest (or highest)
  cost, and returns this information
}
\value{
  The output is a \code{data.frame} with scalar entries
  
  \item{row }{ the row of \code{vacgrid$stops} from the
    \code{"optvac"}-class object corresponding to the lowest
    (or highest) cost }
  \item{col }{ the row of \code{vacgrid$fracs} from the
    \code{"optvac"}-class object corresponding to the lowest
    (or highest) cost }
  \item{frac }{ the actual optimal fraction to vaccinate
    from the \code{vacgrid} table }
  \item{stop }{ the actual optimal (stopping) threshold
    from the \code{vacgrid} table }
  \item{cost }{ the cost associated with the optimal static
    vaccination strategy }
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
\seealso{ \code{\link{optvac}}, \code{\link{plot.optvac}},
  \code{\link{MCepi}} }
\examples{
## for an example of the usage of this function,
## please see the documentation for optvac as
## referenced in the See Also section, above
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
