\name{getcost}
\alias{getcost}
\alias{getvac}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract the Final Costs/Number of Vaccinations of the
  (Managed) Epidemic }
\description{
  Extracts the (summary of the distribution of) final cost(s) or number
  of vaccinations performed of Monte Carlo epidemic under a static
  vaccination strategy as implemented by \code{\link{MCepi}} or a
  single adaptively managed epidemic as implemented by
  \code{\link{manage}} or \code{\link{MCmanage}}
}
\usage{
getcost(obj)
getvac(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ a object of class \code{"MCepi"}, \code{"epiman"},
    or \code{"MCmanage"} }
}
\details{
  \code{getcost} returns distribution information for the
  final (cumulative) costs of the (managed) epidemic,
  and \code{getvac} returns distribution information for
  the cumulative number of vaccinations
  performed.
  
  If \code{x} is if class \code{"epiman"} then the single final
  cost or number of vaccinations is performed is returned.
  
  For further details of this function please see the
  documentation for \code{\link{MCepi}} and \code{\link{manage}}
}
\value{
  The return value is a double-precision scalar indicating the
  median final cost obtained under \code{\link{manage}},
  or a \code{data.frame} giving the median and quantiles calculated
  by \code{\link{MCepi}} or \code{\link{MCmanage}}
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
\seealso{ \code{\link{MCepi}}, \code{\link{manage}},
  \code{\link{MCmanage}} }
\examples{
## for an example of the usage of this function,
## please see the documentation for MCepi as
## referenced in the See Also section, above
}
\keyword{ misc }
