\name{processBAM}
\alias{processBAM}
\alias{processBAMfordups}
\alias{processBAMbylane}
\alias{processBAMSE}
\title{
Functions to process BAM files and compile tables of duplicate fragment sets.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
processBAMfordups(bamfilelist, path, snplist, chromcol = 1, poscol = 3, maxdups = 7, AFcutoff = 0.4, maxdepth = 200, mindepth = 0)

processBAMbylane(bamfilelist, path, snplist, lanelist, chromcol = 1, poscol = 3, maxdups = 7, AFcutoff = 0.4, maxdepth = 200, mindepth = 0)

processBAMSE(bamfilelist, path, snplist, chromcol = 1, poscol = 3, maxdups = 7, AFcutoff = 0.4, maxdepth = 200, mindepth = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamfilelist}{
The (indexed and duplicate-marked) BAM files to process
}
  \item{path}{
Where the BAM files are located
}
  \item{snplist}{
The locations to interrogate
}
\item{lanelist}{
Which lanes of sequencing should be included?
}
  \item{chromcol}{
Which column of snplist has the chromosomes?
}
  \item{poscol}{
Which column of snplist has the position
}
  \item{maxdups}{
How many duplicates in a set so we want to go up to?
}
  \item{AFcutoff}{
What allele fraction (<0.5) do we require to accept a site as being heterozygous?
}
  \item{maxdepth}{
At what maximum depth do we feel that the copy number may be 'wrong' and thus should be excluded
}
  \item{mindepth}{
  At what minimum depth do we feel that the copy number may be 'wrong' and thus should be excluded
}
}
\value{
These functions produce tables of sets of duplicate sequencing fragments catagorized by the allele patterns observed.
}
\references{
Lynch Smith Eldridge and Tavare 2016
}
\author{
Andy Lynch
}

\examples{
## No BAM files are supplied with the package, so we cannot have executable examples

## See the vignette for example usage
}
