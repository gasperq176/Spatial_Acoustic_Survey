\name{scr.sim}
\alias{scr.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Capture Data
}
\description{
Randomly generates a spatial capture-recapture data set, and returns a matrix of encounter rates for each animal at each trap (omega).
}
\usage{
scr.sim(lambda0, sigma, traps, density = 50,
    distr = c("pois", "binom", "nbinom"),
    limits = list(xlim = NULL, ylim = NULL), draw = FALSE,
    binom.n = 1, acoustic = FALSE, lambda_c = NULL,
    toa = FALSE, sigma_toa = 0.002, speed_sound = 330, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda0}{
	Baseline encounter rate (i.e. encounter rate for d = 0)
  }
  \item{sigma}{
    Rate of decay of encounter rates (as distance between trap and animal activity centre increases).
  }
  \item{traps}{
    Matrix of trap location coordinates.
  }
  \item{density}{
	  Point (animal) density in animals per hectare. Default value is 50.
  }
  \item{distr}{
    Distribution of count data (for each trap) to be simulated. Default is Poisson; see 'Details'.
  }
  \item{limits}{
	  List of length 2, giving x and y coordinates ranges. Default value is \code{NULL}; see 'Details'.
  }
  \item{draw}{
    Plot the randomly generated points and traps. Default value is \code{FALSE}.
  }
  \item{binom.n}{
    Specifies the number of trials when generating from the binomial distribution. Default value is 1 (i.e. Bernoulli trial).
  }
  \item{acoustic}{
    Simulate acoustic capture histories. Default value is \code{FALSE}. See 'Details'.
  }
  \item{lambda_c}{
    Mean number of calls made, for acoustic capture data.
  }
  \item{toa}{
    Logical; simulate a time of arrival matrix.
  }
  \item{sigma_toa}{
    Error in recording time. Default value is 0.002.
  }
  \item{speed_sound}{
    Speed of sound. Default is 330(m/s).
  }
  \item{...}{
	  Further \link[graphics:par]{graphical parameters} may also be supplied as arguments.
  }
}
\details{
Total number of points (animal activity centres) is randomly generated from a Poisson distribution: N ~ Pois(DA), where D = \code{dens}, and A = area (automatically calculated from the trap locations; \code{xlim}, \code{ylim} can be specified manually).

Bounding box is calculated from the range given by the coordinates in \code{traps}; points are randomly generated within this area, but not all points will be detected by the traps (depends on \code{sigma}).

Size parameter should be supplied when \code{distr = nbinom}, otherwise default is 2.

Acoustic capture histories contain 1 or more rows per 'call' detected by each trap. The final column of the capture history matrix identifies each animal. Requires \code{lambda_c} (mean number of calls per animal).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum likelihood methods for capture--recapture studies. \emph{Biometrics} \strong{64}, 377--385.
}
\author{
Callum Young
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[secr]{secr}}
\code{\link[scr:scr.fit]{scr::scr.fit}}
}
\examples{
## Simulating with lambda0 = 5, sigma = 15
## Trap locations are stored in data(scr.data)
scr.sim(5, 15, scr.data$traps)

## Acoustic simulation
## Note that for acoustic captures, 0 < lambda0 < 1
scr.sim(1, 10, scr.data$traps, acoustic = TRUE, lambda_c = 3)

%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
