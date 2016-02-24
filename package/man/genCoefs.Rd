\name{genCoefs}
\alias{genCoefs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to generate the multiplicative coeffiecients associated with the probabilities of partitions arising.
}
\description{
Every partition has probability of occurring of the form a*PD^b*FD^c, and this function identifies the  values of a in a recursive manner.
}
\usage{
genCoefs(PL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PL}{
A partitionlist object produced by \code{\link{genPL}}}
}
\details{
Following the methods described in Lynch, Smith et al.
}
\value{
A list of vectors of coefficients, each vector corresponding to a matrix in the partition list.
}
\references{
See the sweave in this package as well as Lynch, Smith et al.
}
\author{
Lynch
}

\seealso{
\code{\link{genPL}}
}
\examples{
## see vignette
}
