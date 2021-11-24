\name{simreads}
\alias{simreads}
\alias{simMTreads}
\title{
Function to simulate sequencing and construct tables of duplicate read-set allele patterns}
\description{
Simulating artificial genomes or mitochondria as required.
}
\usage{
simreads(readlength = 100, depth = 50, pcrdup = 0.05, ISmean = 310, ISsd = 7, realISD = NULL, genlength = 999999, spacing = 1000, lowerlength = 200, upperlength = 400, maxM = NULL)

simMTreads(readlength=100, depthpercopy=25, numcopies, pcrdup=0.05, ISmean=310,ISsd=7, realISD=NULL, genlength=16500, spacing=1000, upperlength=400)
}
\arguments{
  \item{readlength}{
The lengths of the reads (that will be simulated in pairs).
}
  \item{depth}{
The depth of sequencing to be simulated
}
  \item{depthpercopy}{
The depth of sequencing associated with one copy of a genome
}
  \item{numcopies}{
The number of copies of the mitochondrial genome
}
  \item{pcrdup}{
The PCR duplication rate (proportion)
}
  \item{ISmean}{
The insert size mean
}
  \item{ISsd}{
The insert size standard deviation
}
  \item{realISD}{
An insert size distribution from which to sample insert sizes
}
  \item{genlength}{
The genome length
}
  \item{spacing}{
The distance between heterozygous SNPs
}
  \item{lowerlength}{
A lower bound on insert size below which values should be rejected
}
  \item{upperlength}{
An upper bound on insert size below which values should be rejected
}
\item{maxM}{
If specified, and greater than the observed maximum partition size, this is the largest set of duplicate read allele patterns that will be reported
}
}
\value{
These functions return a table of duplicate read-set allele patterns
}
\references{
Lynch Smith Eldridge and Tavare 2016
}
\author{
Andy Lynch
}
\seealso{
\code{\link{processBAMfordups}}
}
\examples{
## See Vignette
}