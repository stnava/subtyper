# subtyper

[![Build Status](https://travis-ci.org/stnava/subtyper.png?branch=master)](https://travis-ci.org/stnava/subtyper)

## Subtype and staging analyses in disease modeling

Common tasks for subtyping include making measurements of data consistency,
*filtering* for good subjects, training the subtyping algorithm, predicting
the subtypes in new data and visualizing results, often over time.

This package expects population-level data frames with longitudinal data.


# Overview of `subtyper`

Subtyper is a package for defining subtypes in patient populations based on data.
A "subtyping method" uses data to classify subjects into groups that represent
different disease etiologies, different progression paths or different sets of
symptoms potentially arising from the same underlying pathology.  Subtyping
therefore tries to identify smaller but distinctive groups within an overall
disorder.

Subtyping methods must be predictive.  That is, they follow a machine learning
style with a "training" phase for the models followed by inference wherein we
apply the fixed models to new data.  Subtyper's design follows this paradigm
which will be seen in the examples below.  While we only use linear regression
here, the basic principles can extend to other modeling techniques including
deep learning.

Subtyping and "staging" are or may be inter-related.  For example, some sub-types
may only be identifiable at a given disease stage.  Other subtypes are fundamentally
different across all disease stages (e.g. genetic forms of AD or PD vs idiopathic forms).
In PD, some patients "jump" between classically defined subtypes (see 2021 Lancet
  review on PD diagnosis).
Ths package does not claim to produce "good" subtypes.  It's just a tool that
allows one to test and visualize different approaches to subtyping whether they
are stage-specific or not.

## Generic processing steps

The brief code block below shows us the core steps in `subtyper`.  We will
comment on each step within its own section.

```{r}
library(subtyper)
mydf = generateSubtyperData( 100 )
rbfnames = names(mydf)[grep("Random",names(mydf))]
mydf = outlierness( mydf, rbfnames )
mydf = highestQualityRepeat( mydf, "Id", "visit", "OL_KNN_SUM")
qdf = trainSubtypeUni( mydf, "cognition", c("C0","C1","C2"), c(0.33,0.66) )
pdf = predictSubtypeUni( mydf, qdf, "Id" )
```

## Example data

The function:

```{r}
mydf = generateSubtyperData( 100 )
```

simply supplies a dataframe that fits our assumptions about how data is organized.

* subject IDs are provided

* each subject has one or more "visits"

* at each visit, there may be more than one repeats of a measurement

* there is some measurement of the quality of each measurement (optional)

While our simulated data is complete, in reality,
many subjects may be missing data.   We can further simulate this  (at random) by
performing:

```{r}
nrows = nrow( mydf )
mixna = sample(1:nrows,10)
mydf[ mixna, "cognition" ] = NA
```

## Outlierness

While we assume that the data comes with some quality measurements (domain specific),
we may also use data-driven approaches to measuring how close each timepoint-repeat
is to related data.  If such a datapoint is "far away" from its nearest neighbors,
then it is likely to be an outlier.  The `outlierness` function implements
several methods for making such estimations which can guide data inspection.

```{r}
rbfnames = names(mydf)[grep("Random",names(mydf))]
mydf = outlierness( mydf, rbfnames )
```

These outlierness measurements are added to the data frame by the above call.
The names are prepended with "OL_".  Which particular measurement -- or set of
measurements -- is best may depend on your dataset.

It is important to note that these methods are likely to be more generalizable
if they are trained and applied to data that is "good" quality.

**Warning:** Some subjects will always be outliers according to these methods.
However, that does not mean they should be rejected.  Decisions about data
rejection should be put off to the last possible moment.  Another option is to
simply covary for the outlierness score --- with the assumption that the score
is not or is only marginally related to the outcomes of interest.

## Filter for quality

Given the above comments, we might want to find the "best" repeat for each
time point in the case when there are multiple repeats.  Here, we use the
`OL_KNN_SUM` score as an outlierness measurement.

```{r}
mydf = highestQualityRepeat( mydf, "Id", "visit", "OL_KNN_SUM")
```


## Covariate adjustment

The function `adjustByCovariates` can be used to train and predict
covariate adjusted scores.  Here, we train the adjustment based on group `G0`
and adjust with respect to `RandomBasisProjection01`.  In practice, these
variables may include age, sex and other nuisance variables.

```{r}
myform = "cognition ~ RandomBasisProjection01 "
mydf = adjustByCovariates(mydf,myform,"DX","G0")
```

## Train the subtyping model

We use a four group model here in order to demonstrate the
difference from the known three group diagnosis.

```{r}
qdf = trainSubtypeUni( mydf, "cognition_adjusted", c("C0","C1","C2","C3"), c(0.25,0.5, 0.75) )
```

## Perform inference with the subtyping model

We define the subtype from the baseline data.  This means that we assume that the
data at baseline is sufficient to confidently identify which grouping to which
the subject should belong.

```{r}
pdf = predictSubtypeUni( mydf, qdf, "Id", "visit", "V0" )
```

## Check subtypes against diagnosis

Data-driven subtypes should overlap --- in our simulated data example --- with diagnosis.

```{r}
knitr::kable( table( pdf$subtype, mydf$DX ), caption='Subtypes vs diagnosis: By design, these categories are very similar.', table.attr = "style='width:60%;'")
```

## Visualize the results

See the difference with diagnosis.

```{r}
summ = plotSubtypeChange( mydf, "Id", "cognition", "DX", "visit", whiskervar='se'  )
```
See the difference with subtype.

```{r}
summ = plotSubtypeChange( pdf, "Id", "cognition", "subtype", "visit", whiskervar='se' )
```

Visualization of the subtype definition against various outcomes is key to
establishing validity.


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
