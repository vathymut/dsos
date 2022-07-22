
# `D-SOS`: Dataset shift with outlier scores

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle)
[![License:
GPL3](https://img.shields.io/badge/License-GPL3-green.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![CRAN](https://www.r-pkg.org/badges/version/dsos)](https://cran.r-project.org/package=dsos)
[![UAI
2022](https://img.shields.io/badge/paper-UAI-yellow)](https://openreview.net/forum?id=S5UG2BLi9xc)
[![Downloads](https://cranlogs.r-pkg.org/badges/dsos)](https://cran.r-project.org/package=dsos)
<!-- badges: end -->

## Overview

`dsos` tests for no adverse shift based on outlier scores. Colloquially,
these tests check whether the new sample is not substantively worse than
the old sample, not if the two are equal as tests of equal distributions
do. `dsos` implements a family of two-sample comparison which assumes
that we have both a training set, the reference distribution, and a test
set.

## Installation

If the package is on [CRAN](https://CRAN.R-project.org), install with

``` r
install.packages("dsos")
```

From GitHub, install with:

``` r
# install.packages("remotes")
remotes::install_github("vathymut/dsos")
```

## Example

Use `dsos` to test for adverse shift on the
[`iris`](https://en.wikipedia.org/wiki/Iris_flower_data_set) dataset.
Here, the outlier scores are from extended isolation forest for
density-based out-of-distribution (OOD) detection:

``` r
library(dsos)
set.seed(12345)
data(iris)
versicolor <- iris[51:100,1:4] # Training sample: Species == 'versicolor'
virginica <- iris[101:150,1:4] # Test sample: Species == 'virginica'
iris_test <- od_pt(x_train = versicolor, x_test = virginica)
plot(iris_test)
```

<img src="man/figures/README-example-1.png" width="100%" />

Among others, `dsos` also implements a method for confidence-based OOD
detection via prediction (resampling) uncertainty.

## Reference

To cite this work and for technical details, please refer to the [arXiv
paper](https://openreview.net/forum?id=S5UG2BLi9xc). Sample Bibtex is
given below:

``` bibtex
@misc{kamulete2021test,
      title={Test for non-negligible adverse shifts}, 
      author={Vathy M. Kamulete},
      year={2021},
      eprint={2107.02990},
      archivePrefix={arXiv},
      primaryClass={stat.ML}
}
```
