\name{libCompNewParam}
\alias{libCompNewParam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
complexity estimate
}
\description{
A function for estimating library complexity. If the complexity proposed is correct, then the value will be zero. 
}
\usage{
libCompNewParam(X, R, C)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
The number of unique molecules in the library
}
  \item{R}{
Duplication Rate
}
  \item{C}{
Number of unique molecules sequenced
}
}
\details{
See Sweave for more details
}
\value{
A measure of discordence. This will be zero if the value of X is correct.
}
\references{
See the sweave in this package as well as Lynch, Smith et al.
}
\author{
Lynch
}
\examples{
## see vignette
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ complexity }
