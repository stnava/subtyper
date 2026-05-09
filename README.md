# subtyper

[![CircleCI](https://circleci.com/gh/stnava/subtyper/tree/master.svg?style=svg)](https://circleci.com/gh/stnava/subtyper/tree/master)

documentation page [here](https://stnava.github.io/subtyper/)

## Subtype and staging analyses in disease modeling

Common tasks for subtyping include making measurements of data consistency,
filtering for subjects with good/high quality data, training the subtyping 
algorithm, predicting the subtypes in new data and visualizing results, 
often over time.  Many helper functions of data exploration, manuscript 
preparation, model checking, covariate or site/confound control etc.

This package expects population-level data frames with longitudinal data.


## Installing

The pre-release version of the package can be pulled from GitHub using [pak](https://pak.r-lib.org/) (recommended) or [remotes](https://remotes.r-lib.org/):

```r
    # Using pak (handles all dependencies automatically)
    pak::pkg_install("stnava/subtyper")

    # Or using remotes
    remotes::install_github("stnava/subtyper", build_vignettes=TRUE)
```

**Note on Archived Dependencies:** This package depends on `imbalance` and `DDoutlier`, which have been archived on CRAN. They are automatically installed from their GitHub mirrors when using `pak` or `remotes` because they are listed in the `Remotes:` field of the `DESCRIPTION` file.

See the file `.circleci/config.yml` for hints about system-level dependencies.


## OSX specific notes: Running R in a Clean Environment

### Why Use This?

Sometimes, issues arise in R installations due to user-specific environment variables (such as `LD_LIBRARY_PATH`, `DYLD_LIBRARY_PATH`, `R_LIBS`, or `PKG_CONFIG_PATH`). These can conflict with R's ability to find the correct system libraries or load shared objects (`.so` / `.dylib` files) during package compilation or loading.

To avoid these issues, you can launch R in a **minimal, clean environment**, preserving only the necessary `HOME` and `PATH` variables.

### Command

```bash
env -i HOME=$HOME PATH="/opt/homebrew/bin:/usr/bin:/bin" /opt/homebrew/bin/R
```

## For developers

The repository includes a Makefile to facilitate some common tasks.

### Running tests

`$ make test`. Requires the [testthat](http://testthat.r-lib.org/) package. You can also specify a specific test file or files to run by adding a "file=" argument, like `$ make test file=logging`. `testthat::test_package()` will do a regular-expression pattern match within the file names (ignoring the `test-` prefix and the `.R` file extension).

### Updating documentation

`$ make doc`. Requires the [roxygen2](https://github.com/klutometis/roxygen) package.
