# subtyper

[![CircleCI](https://circleci.com/gh/stnava/subtyper/tree/master.svg?style=svg)](https://circleci.com/gh/stnava/subtyper/tree/master)

documentation page [here](https://stnava.github.io/subtyper/)

## Subtype and staging analyses in disease modeling

Common tasks for subtyping include making measurements of data consistency,
filtering for subjects with good/high quality data, training the
subtyping algorithm, predicting the subtypes in new data and
visualizing results, often over time.

This package expects population-level data frames with longitudinal data.


## Installing

The pre-release version of the package can be pulled from GitHub using the [devtools](https://github.com/r-lib/devtools) package:

```r
    # install.packages("devtools")
    devtools::install_github("stnava/subtyper", build_vignettes=TRUE)
```

See the file `.circleci/config.yml` for hints about installing the many dependencies.  A few other hints below.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("Biobase", quietly = TRUE))
  BiocManager::install("Biobase")
install.packages("NMF")
```

```r
sudo apt install libgmp-dev
# then install ClusterR
```

`pgenlibr`

```
# brew install libsm
gh repo clone chrchang/plink-ng
cd plink-ng/2.0/pgenlibr/
# edit src/Makevars to include  -I/System/Volumes/Data/Users/stnava/code/extern/plink-ng/2.0/simde
R CMD INSTALL . # FIXME not on arm64 yet
```

```R
install.packages(c("ClusterR", "effectsize", "Evacluster", "flexclust", "ggthemes", "ggstatsplot", "ggbeeswarm", "wesanderson", "Hmisc", "DDoutlier", "NMF", "pheatmap", "gprofiler2", "pgenlibr","gaston"))
```

## For developers

The repository includes a Makefile to facilitate some common tasks.

### Running tests

`$ make test`. Requires the [testthat](http://testthat.r-lib.org/) package. You can also specify a specific test file or files to run by adding a "file=" argument, like `$ make test file=logging`. `testthat::test_package()` will do a regular-expression pattern match within the file names (ignoring the `test-` prefix and the `.R` file extension).

### Updating documentation

`$ make doc`. Requires the [roxygen2](https://github.com/klutometis/roxygen) package.
