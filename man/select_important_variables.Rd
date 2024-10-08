% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subtyper.R
\name{select_important_variables}
\alias{select_important_variables}
\title{Select important variables based on stabilized correlations}
\usage{
select_important_variables(
  data,
  cols,
  threshold = 0.5,
  epsilon = 0.001,
  max_iter = 0,
  use_glasso = 0.1,
  least = FALSE,
  return_raw_importance = FALSE
)
}
\arguments{
\item{data}{A data frame containing the variables.}

\item{cols}{A character vector specifying which columns (variables) to consider.}

\item{threshold}{Proportion of variables to select based on importance (default: 0.2).}

\item{epsilon}{Convergence threshold for Sinkhorn iterations (default: 1e-3).}

\item{max_iter}{Maximum number of Sinkhorn iterations (default: 100).}

\item{use_glasso}{set a scalar greater than zero to use graphical LASSO; this parameter relates to sparseness levels}

\item{least}{boolean take least explanatory (most uncorrelated) variables}

\item{return_raw_importance}{boolean}
}
\value{
A character vector of selected variable names.
}
\description{
This function selects important variables from a dataset based on stabilized correlations
computed using the Sinkhorn method.
}
\examples{
set.seed(123)
data <- data.frame(
  x1 = rnorm(100),
  x2 = rnorm(100),
  x3 = rnorm(100),
  x4 = rnorm(100)
)
selected_vars <- select_important_variables(data, 
 c("x1", "x2", "x3", "x4"), threshold = 0.3)
print(selected_vars)
}
