\name{propllik}
\alias{propllik}
\title{
Function to estimate the likelihood of a given proportion of duplicate reads being fragmentation duplicates
}
\description{
This function takes in a set of observations of allele patterns, from sites presumed to be diploid heterozygous, and calculates the log-likelihood of the F_D parameter through calls to \code{genprobs}. It is currently hard-coded to deal with sets of up to 7 duplicate fragments.
}
\usage{
propllik(x, FD)
}
\arguments{
  \item{x}{
A vector describing the nuber of citings of each allele pattern, in the order: AA, AB, AAA, AAB, AAAA, AAAB, AABB, AAAAA, AAAAB, AAABB, AAAAAA, AAAAAB, AAAABB, AAABBB, AAAAAAA, AAAAAAB, AAAAABB, AAAABBB.
}
  \item{FD}{
The value of the F_D parameter for which we desire to calculate the likelihood
}
}
\details{
Following the methods described in Lynch, Smith et al.
}
\value{
The log-likelihood of the parameter in question.
}
\references{
See the sweave in this package as well as Lynch, Smith et al.
}
\author{
Lynch
}
\seealso{
\code{\link{genprobs}}
}
\examples{
## see vignette
}
