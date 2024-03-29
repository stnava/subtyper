% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subtyper.R
\name{adjustByCovariates}
\alias{adjustByCovariates}
\title{Covariate adjustment}
\usage{
adjustByCovariates(
  mxdfin,
  adjustmentFormula,
  groupVariable,
  group,
  allowMissing = FALSE
)
}
\arguments{
\item{mxdfin}{Input data frame with repeated measurements and a grouped time variable}

\item{adjustmentFormula}{string defining a valid formula to be used in \code{lm}.}

\item{groupVariable}{names of the column that defines the group to use for training.}

\item{group}{string defining a subgroup on which to train}

\item{allowMissing}{boolean}
}
\value{
data frame with adjusted measurement variable as defined in formula
}
\description{
Adjust a training vector value by nuisance variables eg field strength etc.
One may want to use a specific sub-group for this, e.g. controls.
}
\examples{
mydf = generateSubtyperData( 100 )
myform = "cognition ~ RandomBasisProjection01 "
mydf = adjustByCovariates(mydf,myform,"DX","G0")
}
\author{
Avants BB
}
