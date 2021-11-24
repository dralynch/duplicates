\name{libCompNewParam}
\alias{libCompNewParam}
\title{
Function to confirm library complexity}
\description{
Function to confirm library complexity following approach of Picard
}
\usage{
libCompNewParam(X, R, N)
}
\arguments{
  \item{X}{
Complexity
}
  \item{R}{
Duplicate rate
}
  \item{N}{
Number of molecules sequenced
}
}
\value{
0 if complexity estimate is 'correct', higher otherwise.
}
\references{
Lynch Smith Eldridge and Tavare 2016
}
\author{
Andy Lynch
}
\examples{
## See vignette
}