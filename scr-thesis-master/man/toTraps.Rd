\name{toTraps}
\alias{toTraps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Format matrix of traps co-ordinates
}
\description{
Convert a matrix of trap co-ordinates into a format expected by \code{read.traps()}.
}
\usage{
toTraps(traps)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{traps}{
	  Matrix of x- and y- co-ordinates.
  }
}
\details{
Detector IDs are automatically given as \code{1:nrow(traps)}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A data frame is returned, with detector IDs, and co-ordinates labelled "x" and "y".

This data frame can be read by \code{secr::read.traps()}.
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
\code{\link[secr:read.traps]{secr::read.traps}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Load the data exported with scr
## Then convert the example capture history

load(scr.data)

read.traps(data = toTraps(scr.data$traps), detector = "count")

%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
