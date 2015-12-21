coala 0.3.0
===========

## Major improvements
* Support for more selection models, including ones for local adaptation (#137).
* Add `as.segsites.GENOME` function that converts genetic data imported with
  the package PopGenome to coala's format (#139).

## Small Changes
* Adds `feat_ignore_singletons`, which is a feature that makes coala ignore 
  singletons when calculating the summary statistics (#138).
* Use `ms` from package `phyclust` instead of requiring that the binary is
  installed on the system (#140).
* Ensure that msms uses only one CPU core (#142).



coala 0.2.2
===========

* Fixes the broken nucleotide diversity and Tajima's D summary statistics
  (#133).
* Adds support for calculating joint frequency spectra for more than two 
  populations (#132).
  


coala 0.2.1
===========

* Fixes a test that failed on R 3.1.x due to a bug in the tests code (#127).
* Fixes version requirement for `testthat` (#127).
* Adds the `calc_sumstats_from_data` function for calculating summary statistics 
  from biological data (#124).
* Exports the functions related to segregating sites (#122).



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
