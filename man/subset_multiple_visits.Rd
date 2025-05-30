% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subtyper.R
\name{subset_multiple_visits}
\alias{subset_multiple_visits}
\title{Subset Dataframe for Subjects with Multiple Unique Visits}
\usage{
subset_multiple_visits(data, subject_identifier, visit_identifier)
}
\arguments{
\item{data}{A dataframe containing the data to be subset.}

\item{subject_identifier}{A string specifying the column name used to identify subjects.}

\item{visit_identifier}{A string specifying the column name used to identify visits.}
}
\value{
A dataframe containing only subjects with more than one unique visit.
}
\description{
This function subsets a dataframe to return only subjects who have more than one unique visit.
It takes in the column names for the subject identifier and the visit identifier, then filters
the dataframe to include only rows corresponding to subjects with multiple distinct visits.
}
\examples{
# Example dataframe
df <- data.frame(
  subjectID = c(1, 1, 2, 2, 3, 3, 4),
  visitID = c(1, 2, 1, 1, 1, 3, 1)
)

# Subset for subjects with multiple unique visits
result <- subset_multiple_visits(df, "subjectID", "visitID")
print(result)

}
