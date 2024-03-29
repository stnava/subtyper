% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subtyper.R
\name{mapLRAverageVar}
\alias{mapLRAverageVar}
\title{Convert left/right variables to an average measurement}
\usage{
mapLRAverageVar(
  mydataframe,
  leftvar,
  leftname = "left",
  rightname = "right",
  replacer = "LRAVG"
)
}
\arguments{
\item{mydataframe}{dataframe containing relevant variables}

\item{leftvar}{left side variable names ie the full names of the variables to average}

\item{leftname}{the variable substring indicating left side}

\item{rightname}{the variable substring indicating right side}

\item{replacer}{string to replace left with in column names of output}
}
\value{
fixed x
}
\description{
Convert left/right variables to an average measurement
}
\author{
Avants BB
}
