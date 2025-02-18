# MultiTraits 0.0.2 (2024-12-21)

## Bug fixes

* Fixed bug where discrete scales could not map aesthetics only consisting of
  `NA`s (#5623)
* Fixed spurious warnings from `sec_axis()` with `breaks = NULL` (#5713).
* Patterns and gradients are now also enabled in `geom_sf()` 
  (@teunbrand, #5716).
* The default behaviour of `resolution()` has been reverted to pre-3.5.0 
  behaviour. Whether mapped discrete vectors should be treated as having 
  resolution of 1 is controlled by the new `discrete` argument.
* Fixed bug in `guide_bins()` and `guide_coloursteps()` where discrete breaks,
  such as the levels produced by `cut()`, were ordered incorrectly 
  (@teunbrand, #5757).
  
## Improvements

## Bug Fixes

* Fixed calculation issues in CSR module

* Resolved NaN errors in TN module

## Features

* Optimized visualization effects across all modules (CSR, LHS, NPT, TN)

## Data

* Removed WH dataset

* Updated examples across package modules
