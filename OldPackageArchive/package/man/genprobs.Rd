\name{genprobs}
\alias{genprobs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to work out the probabilities of different allele patterns given the fragmentation duplicate proportion.
}
\description{
This function works out the probabilities of seeing specific allele patterns given assignments of read-pairs to original molecules, and the probabilities of seeing patterns of assignments of read-pairs to molecules given the fragmentation duplicate proportion in order to generate the probability of the allele patterns given the fragmentation duplicate proportion.
}
\usage{
genprobs(N, FD, fullout = F, quickcoef = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
The number of duplicate read-pairs in the set being conisdered.
}
  \item{FD}{
The proportion of duplicate read-pairs that are fragmentation duplicates
}
  \item{fullout}{
If true, then the output consists of the partitions, the parameters of the partitions, the probabilities of the partitions, and the probabilities of the allele patterns
}
\item{quickcoef}{
If False, then an object called CoefList is assumed to exist, containing the binary coefficients for the possible partitions of read-pairs into originating molecules. If True, then a quick (but harder to explain) calculation is performed.
}
}
\details{
Following the methods described in Lynch, Smith et al.
}
\value{
a vector of probabilities for the potential allele frequencies associated with \code{N} read-pairs.
}
\references{
Following the methods described in Lynch, Smith et al.
}
\author{
Lynch
}
\seealso{
\code{\link{genCoefs}}, 
}
\examples{
## see vignette
}
