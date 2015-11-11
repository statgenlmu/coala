coala 0.2.0
===========

## Major improvements
* Adds support for distribution independent repetitions on multiple CPU 
  cores (#116).
* Improves support for polyploid models. The `ploidy` parameter is now
  provided in the `coal_model` instead of in `feat_unphased` (#115).
* Adds the MCMF summary statistic (#94).
* Adds support for the omega statistic using OmegaPlus (#109).

## Small Changes
* Adds option to calculate iHS in `sumstat_ihh()` and made the statistic return
  a `data.frame` instead of a list. 
* Adds optional support for calculating the JSFS per locus instead of 
  globally (#112).
* Adds optional in-place transformation of summary statistics (#110).
* Adds support for simulating a fixed number of mutations with ms and 
  msms (#19).
* Writes seq-gen output into memory instead of in files before it is 
  parsed (#99).
* Adds optional support for parameter zero inflation for a deterministic
  fraction of loci instead of a random number. Can be used by setting
  `random = FALSE` in `par_zero_inflation` (#97).
* `get_outgroup` now returns `NA` if the model has no outgroup rather than 
  throwing an error.

## Bug Fixes
* Fixes the simulation of sizes changes in one populations models with msms (#105).
* Remove broken implementation of the nSL statistic (`sumstat_nsl()`)
* Fixes site frequency calculation when an outgroup is present (#96).
* Fixes multiple errors that occurred in edge cases when calculating ihh (#98).



coala 0.1.1
===========

* Fixes a memory corruption that occurred only in tests (#90).
* Updates `README.md`.
* Corrects various typos.



coala 0.1.0
===========

* Initial release version
* Thanks to Ann Kathrin Huylmans for suggesting the name 'coala' and
  to Soumya Ranganathan for proofreading.
