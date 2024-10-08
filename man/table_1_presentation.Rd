% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subtyper.R
\name{table_1_presentation}
\alias{table_1_presentation}
\title{Format a table for publication using LaTeX or html}
\usage{
table_1_presentation(
  df,
  caption = "",
  label = "",
  format = "latex",
  html_font = "Arial"
)
}
\arguments{
\item{df}{Data frame containing the table data.}

\item{caption}{Caption for the table.}

\item{label}{Label for the table.}

\item{html_font}{character eg Arial}
}
\value{
A LaTeX-formatted table.
}
\description{
Format a table for publication using LaTeX or html
}
\examples{
\dontrun{
table_1_presentation(df, caption = "Demographic and Clinical Characteristics",
                         label = "table:demographics")
}
}
