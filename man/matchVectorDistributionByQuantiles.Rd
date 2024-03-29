% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subtyper.R
\name{matchVectorDistributionByQuantiles}
\alias{matchVectorDistributionByQuantiles}
\title{Match a pair of vector distributions based on quantiles}
\usage{
matchVectorDistributionByQuantiles(
  vecToBeTransformed,
  vecReference,
  quantiles = 1:9/10,
  polynomialOrder = 1,
  truncate = TRUE
)
}
\arguments{
\item{vecToBeTransformed}{input vector to be transformed to match the reference}

\item{vecReference}{reference vector}

\item{quantiles}{a vector of quantile points to match}

\item{polynomialOrder}{integer greater than or equal to one}

\item{truncate}{boolean}
}
\value{
the transformed vector
}
\description{
Match a pair of vector distributions based on quantiles
}
\examples{
mydf = generateSubtyperData( 100 )
rbfnames = names(mydf)[grep("Random",names(mydf))]
newvec = matchVectorDistributionByQuantiles( mydf[,rbfnames[1]], mydf[,rbfnames[1]] )
}
\author{
Avants BB
}
