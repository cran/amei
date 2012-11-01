\name{print.Rd}
\alias{print.MCepi}
\alias{print.optvac}
\alias{print.epiman}
\alias{print.summary.MCepi}
\alias{print.summary.optvac}
\alias{print.summary.epiman}
\alias{summary.MCepi}
\alias{summary.optvac}
\alias{summary.epiman}
\title{ Summarizing and printing amei output }
\description{
  Summarizing, printing, and plotting the contents of a
  the following objects  \code{"MCepi"}, \code{"optvac"},
  \code{"epiman"}.  In the current version the printing
  and summary commands are very similar
}
\usage{
\method{print}{MCepi}(x, ...)
\method{print}{optvac}(x, ...)
\method{print}{epiman}(x, ...)
\method{summary}{MCepi}(object, ...)
\method{summary}{optvac}(object, ...)
\method{summary}{epiman}(object, ...)
\method{print}{summary.MCepi}(x, ...)
\method{print}{summary.optvac}(x, ...)
\method{print}{summary.epiman}(x, ...)
}
\arguments{
  \item{object}{ an object of class \code{"MCepi"}, \code{"optvac"}, or
  \code{"epiman"} that must be named
    \code{object} for the generic methods \code{\link{summary}} }
  \item{x}{ an object of class \code{"MCepi"}, \code{"optvac"}, or
    \code{"epiman"} that must be named \code{x} for generic printing
    and plotting via  \code{\link{print}} }
  \item{...}{ passed to the generic \code{\link{print}} or
    \code{\link{summary}} commands}
}
\details{
  The printing and summaries provided by these functions are
  essentially identical except that \code{summary} returns a \code{list}
  that allows the information to be extracted for external use
  (in code) whereas the printing adds some additional text for
  human consumption including information about the call

  Other ways of extracting information contained in these objects
  include the functions \code{\link{getvac}}, \code{\link{getcost}},
  and \code{\link{getpolicy}}.  The plotting functions
  \code{\link{plot.MCepi}}, \code{\link{plot.optvac}},
  and \code{\link{plot.epiman}} are also helpful
  
  The \code{list}(s) returned by the \code{summary} command are
  detailed in the value section below
}
\value{
  The \code{\link{summary}} commands documented here return
  a list containing (a subset of) the items below.
  The other functions do not return values.

  \item{obj}{ a copy of the input object }
  \item{final}{ in the case of \code{"MCepi"}-class object this
    is a \code{data.frame} containing a summary of the distribution
    of the total number of vaccinations and the final cost; for
    \code{"epiman"}-class objects the distributional information
    is replaced by a scalar }
  \item{params}{ in the case of \code{"epiman"}-class objects
    this \code{data.frame} containing a summary of the distributional
    information of SIR model parameters obtained during the final
    time step of the epidemic }
  \item{best}{ the best policy contained in the \code{"optvac"}-class
    object }
  \item{worst}{ the worst policy contained in the \code{"optvac"}-class
    object }
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

\seealso{ \code{\link{MCepi}}, \code{\link{optvac}},
  \code{\link{manage}}, \code{\link{MCmanage}},
  \code{\link{getvac}}, \code{\link{getcost}}, \code{\link{getpolicy}},
  \code{\link{plot.MCepi}}, \code{\link{plot.optvac}},
  \code{\link{plot.epiman}}}
\keyword{ methods }
