# subtyper

[![Build Status](https://travis-ci.org/stnava/subtyper.png?branch=master)](https://travis-ci.org/stnava/subtyper)

documentation page [here](https://stnava.github.io/subtyper/)

## Subtype and staging analyses in disease modeling

Common tasks for subtyping include making measurements of data consistency,
*filtering* for good subjects, training the subtyping algorithm, predicting
the subtypes in new data and visualizing results, often over time.

This package expects population-level data frames with longitudinal data.


## Installing

The pre-release version of the package can be pulled from GitHub using the [devtools](https://github.com/r-lib/devtools) package:

    # install.packages("devtools")
    devtools::install_github("stnava/subtyper", build_vignettes=TRUE)

## For developers

The repository includes a Makefile to facilitate some common tasks.

### Running tests

`$ make test`. Requires the [testthat](http://testthat.r-lib.org/) package. You can also specify a specific test file or files to run by adding a "file=" argument, like `$ make test file=logging`. `testthat::test_package()` will do a regular-expression pattern match within the file names (ignoring the `test-` prefix and the `.R` file extension).

### Updating documentation

`$ make doc`. Requires the [roxygen2](https://github.com/klutometis/roxygen) package.
