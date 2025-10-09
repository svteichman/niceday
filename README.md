
<!-- README.md is generated from README.Rmd. Please edit that file -->

# niceday

<!-- badges: start -->

[![R-CMD-check](https://github.com/statdivlab/niceday/workflows/R-CMD-check/badge.svg)](https://github.com/statdivlab/niceday/actions)
<!-- badges: end -->

`niceday` is an `R` package for estimating fold-differences in
multivariate outcomes across binary treatments and in the presence of
other covariates. It uses a nonparametric approach, which means that
under a given set of identifiability assumptions, no assumptions about
the distribution of outcomes or the relationship between covariates and
outcomes need to be made.

## Installation

To install and load `niceday`, use the following code.

``` r
# install.packages("remotes")
remotes::install_github("statdivlab/niceday")
library(niceday)
```

## Use

Consider a count matrix `W`, metadata `data`, the treatment `A`, and
covariates to adjust for `X`.

``` r
library(niceday)
#> Loading required package: gam
#> Loading required package: splines
#> Loading required package: foreach
#> Loaded gam 1.22-6
#> Loading required package: MASS
#> Loading required package: glmnet
#> Loading required package: Matrix
#> Loaded glmnet 4.1-10
#> Loading required package: hal9001
#> Loading required package: Rcpp
#> hal9001 v0.4.6: The Scalable Highly Adaptive Lasso
#> note: fit_hal defaults have changed. See ?fit_hal for details
#> Loading required package: ranger
#> Loading required package: xgboost
#> Loading required package: withr
#> Loading required package: SuperLearner
#> Loading required package: nnls
#> Super Learner
#> Version: 2.0-29
#> Package created on 2024-02-06

data(EcoZUR_meta)
data(EcoZUR_count)
ndFit(W = EcoZUR_count[, 1:50], # consider only the first 50 taxa to run quickly
      data = EcoZUR_meta,
      A = ~ Diarrhea,
      X = ~ sex + age_months,
      num_crossval_folds = 2, # use more folds in practice
      num_crossfit_folds = 2, # for cross validation and cross fitting
      sl.lib.pi = c("SL.mean"), # choosing single learner for the example to run quickly,
      sl.lib.m = c("SL.mean"))  # in practice would use other options as well
#> 
#> Call:
#> ndFit(W = EcoZUR_count[, 1:50], A = ~Diarrhea, X = ~sex + age_months, 
#>     data = EcoZUR_meta, num_crossval_folds = 2, num_crossfit_folds = 2, 
#>     sl.lib.pi = c("SL.mean"), sl.lib.m = c("SL.mean"))
#> 
#> 
#> Coefficient estimates with the largest magnitudes:
#>    category        est        se  lower_marg  upper_marg  lower_sim  upper_sim
#> 13       13 -3.2741166 0.6118426 -4.47330614 -2.07492704 -5.2715008 -1.2767324
#> 3         3 -2.6157685 0.7393277 -4.06482416 -1.16671283 -5.0293326 -0.2022044
#> 22       22  2.1413861 0.9688809  0.24241444  4.04035780 -1.0215634  5.3043356
#> 39       39  2.0849190 0.7705275  0.57471289  3.59512520 -0.4304980  4.6003361
#> 30       30 -1.8174642 0.5879262 -2.96977845 -0.66515004 -3.7367723  0.1018438
#> 14       14  1.7374620 0.4585679  0.83868546  2.63623844  0.2404493  3.2344746
#> 21       21  1.6905260 0.5946739  0.52498671  2.85606538 -0.2508099  3.6318620
#> 20       20  1.5157064 0.5558763  0.42620878  2.60520394 -0.2989736  3.3303863
#> 33       33 -1.4345277 0.5776259 -2.56665373 -0.30240164 -3.3202101  0.4511547
#> 32       32 -1.3580239 0.6794349 -2.68969173 -0.02635606 -3.5760655  0.8600177
#> 15       15 -1.3004635 0.6041093 -2.48449595 -0.11643105 -3.2726018  0.6716748
#> 31       31  1.2779538 0.6321374  0.03898736  2.51692033 -0.7856833  3.3415909
#> 45       45 -1.2336750 0.4405577 -2.09715223 -0.37019774 -2.6718928  0.2045428
#> 27       27  1.1777250 0.6579759 -0.11188403  2.46733410 -0.9702630  3.3257130
#> 42       42 -1.0418668 0.5415186 -2.10322383  0.01949026 -2.8096755  0.7259420
#> 12       12  1.0119062 0.8021781 -0.56033395  2.58414645 -1.6068355  3.6306480
#> 25       25  1.0117772 0.3666362  0.29318338  1.73037094 -0.1851211  2.2086754
#> 23       23  1.0073781 0.5866781 -0.14248979  2.15724602 -0.9078554  2.9226116
#> 8         8 -0.8832539 0.4902861 -1.84419704  0.07768917 -2.4838120  0.7173042
#> 4         4 -0.8691163 0.3683387 -1.59104691 -0.14718560 -2.0715724  0.3333399
#>    type
#> 13 TMLE
#> 3  TMLE
#> 22 TMLE
#> 39 TMLE
#> 30 TMLE
#> 14 TMLE
#> 21 TMLE
#> 20 TMLE
#> 33 TMLE
#> 32 TMLE
#> 15 TMLE
#> 31 TMLE
#> 45 TMLE
#> 27 TMLE
#> 42 TMLE
#> 12 TMLE
#> 25 TMLE
#> 23 TMLE
#> 8  TMLE
#> 4  TMLE
#> To obtain the entire coefficient table, use the command `ndFit_object[[1]]$coef` for the Psi1 or Psi2 parameter and `ndFit_object[[2]]$coef` for the Psi1g or Psi2g parameter.
```

## Citation

If you use `niceday` for your analysis, please cite our manuscript.

Grant Hopkins, Sarah Teichman, Ellen Graham, and Amy Willis.
“Nonparametric Identification and Estimation of Ratios of Multi-Category
Means under Preferential Sampling.”

## Bug reports and feature requests

If you identify a bug in `niceday` or have a feature request, please
post an issue [here](https://github.com/statdivlab/niceday/issues).

## Nomenclature

`niceday` stands for Nonparametrically Identified (de-Confounded)
Estimands for Difference Abundance, Yippee!
