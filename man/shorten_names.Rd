% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subtyper.R
\name{shorten_names}
\alias{shorten_names}
\title{Shorten Names}
\usage{
shorten_names(names, max_length = 20, custom_replacements = NULL)
}
\arguments{
\item{names}{A vector of names to shorten.}

\item{max_length}{The maximum length of the shortened names. Defaults to 20.}

\item{custom_replacements}{A list of custom replacements to apply. Defaults to NULL.}
}
\value{
A vector of shortened names.
}
\description{
Shorten a vector of names by removing common substrings, applying custom replacements,
replacing "__" with "\emph{", then "}" with ".", removing duplicates, and truncating to a maximum length.
}
\examples{
names <- c("npsy_BDI_Total", "npsy_BSI.18_TotalRaw", "npsy_CVLTShortDelayFreeRecall_Raw")
shortened_names <- shorten_names(names)
print(shortened_names)
}
