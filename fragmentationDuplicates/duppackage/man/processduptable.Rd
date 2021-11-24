\name{processduptable}
\alias{processduptable}
\alias{genprobsmat}
\title{
Function to process tables of duplicate read-set allele patterns
}
\description{
Takes a table of duplicate read-set allele patterns and returns for each row the observed duplicate rate (percentage), the proportion of duplicates that are due to fragmentation, the corrected PCR duplicate rate (percentage), and the duplicate rate due to fragmentation (percentage)
}
\usage{
processduptable(HST, PL = NULL, CL = NULL, reso = 1001)
genprobsmat(N,FD,PL,CL,reso=1001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{HST}{
A table of duplicate read-set allele patterns such as that produced by processBAMfordups
}
  \item{PL}{
A pre-compiled partition list as produced by genPL() (will be generated if not provided)}
  \item{CL}{A pre-compiled list of coefficients as produced by genCoefs() (will be generated if not provided)}
  \item{reso}{
The resolution to which we wish to know the fragmentation duplicate proportion (1/(reso-1) is achieved)
}
  \item{N}{
  The size of the duplicate read-set being considered
  }
  \item{FD}{
  The Fragmentation duplicate proportion being considered.
  }
}
\value{
A table with four columns: The observed duplication rate (percentage), the proportion of duplicates due to fragmentation, the estimated (corrected) PCR duplication rate (percentage), the estimated duplication rate (percentage) due to coincidental fragmentation.
}
\references{
Lynch Smith Eldridge and Tavare 2016
}
\author{
Andy Lynch
}
\seealso{
processBAMfordups
}
\examples{
## See vignette
}
