\name{toCapthist}
\alias{toCapthist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert matrix into capthist format
}
\description{
Convert a matrix returned by \code{scr.sim()} into a format expected by \code{make.capthist()}.
}
\usage{
toCapthist(captures)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{captures}{
	Matrix returned by \code{scr.sim()}.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A data frame is returned with \code{session}, \code{ID}, \code{occasion}, and \code{trap} columns. The data frame is formatted such that it may be read by \code{\link[secr:make.capthist]{secr::make.capthist()}}. Each row corresponds to a capture, and describes the session and trapping occasion that it was detected during.
}
\references{
Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum likelihood methods for capture--recapture studies. \emph{Biometrics} \strong{64}, 377--385.

Efford, M. G., Borchers D. L. and Byrom, A. E. (2009) Density estimation by spatially explicit capture-recapture: likelihood-based methods. In: D. L. Thomson, E. G. Cooch and M. J. Conroy (eds) \emph{Modeling Demographic Processes in Marked Populations}. Springer, New York. Pp. 255--269.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Load the data exported with scr
## Then convert the example capture history

load(scr.data)

toCapthist(scr.data$capt)

%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
