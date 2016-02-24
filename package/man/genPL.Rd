\name{genPL}
\alias{genPL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to generate a list of partitions.
}
\description{
This function recursively generates a list of the possible partitions of a set of N objects for N=1 up to N=LPL}
\usage{
genPL(LPL = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LPL}{
The upper limit on the number of objects to consider. Don't increase it unless you have read the sweave!
}
}

\value{
A list of matrices, the rows of the Nth matrix providing the partitions N read-pairs
}
\references{
See the sweave in this package as well as Lynch, Smith et al.
}
\author{
Lynch
}
\seealso{
\code{\link{genCoefs}}
}
\examples{
## see vignette
}