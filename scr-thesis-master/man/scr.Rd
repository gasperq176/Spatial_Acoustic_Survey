\name{scr}
\alias{scr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spatial Capture-Recapture Models
}
\description{
Functions to simulate and fit models to spatial capture-recapture data that has been sampled using an array of passive detectors. Models estimate animal population density, and other such parameters related to SCR.
}
\details{
\tabular{ll}{
Package: \tab scr\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2018-02-15\cr
License: \tab MIT\cr
LazyLoad: \tab yes\cr
}

%Spatial capture-recapture refers to the family of methods used to study populations of animals distributed throughout some space. The data are given as capture histories,

}
\references{
Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum likelihood methods for capture--recapture studies. \emph{Biometrics} \strong{64}, 377--385.

Efford, M. G., Borchers D. L. and Byrom, A. E. (2009) Density estimation by spatially explicit capture-recapture: likelihood-based methods. In: D. L. Thomson, E. G. Cooch and M. J. Conroy (eds) \emph{Modeling Demographic Processes in Marked Populations}. Springer, New York. Pp. 255--269.
}
\author{
Callum Young \email{cyou878@aucklanduni.ac.nz}
}
\acknowledgements{
This package absolutely wouldn't exist (or work) without Ben Stevenson. Without Murray Efford's work on SCR and the \code{secr} package,
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[secr]{secr}}
\code{\link[scr:scr.sim]{scr.sim}}
\code{\link[scr:scr.fit]{scr::scr.fit}}
\code{\link[scr:scr.data]{scr::scr.data}}
}
\examples{
\dontrun{
## Simulating with lambda0 = 5, sigma = 15
## Trap locations are stored in data(scr.data)
scr.sim(5, 15, scr.data$traps)

## Acoustic simulation
## Note that for acoustic captures, 0 < lambda0 < 1
scr.sim(1, 10, scr.data$traps, acoustic = TRUE, lambda_c = 3)

## Model fit
scr.fit(scr.data$bin.capt, scr.data$traps, scr.data$mask)

%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ scr}% use one of  RShowDoc("KEYWORDS")
