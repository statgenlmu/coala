coalsimr
========

[![Build Status](https://travis-ci.org/paulstaab/coalsimr.png?branch=master)](https://travis-ci.org/paulstaab/jaatha) 
[![Coverage Status](https://coveralls.io/repos/paulstaab/coalsimr/badge.svg?branch=master)](https://coveralls.io/r/paulstaab/coalsimr)

This is a framework for calling coalescent simulators from within R. It allows to 
specify a model using a ggplot2-like syntax and supports multiple simulators. It can
also incorporate the program _seq-gen_ to simulate finite site mutation models.

This package is currently in early development. An initial CRAN release will occur when the
milestone `0.1` is reached.

The development version can be installed via

```R
devtools::install_github('paulstaab/coalsimr')
```

Use it at your own risk. Please mind that the documentation is currently incomplete and 
partly outdated, and that the user interface will change frequently. 


Supported Simulators
--------------------
The package supports the coalescent simulators

* _ms_,
* _scrm_ and
* _msms_.

All simulators can be combined with _seq-gen_ to simulate finite sites mutation models.
The programs _msms_ and _seq-gen_ must be installed manually on the system, while the
R versions of _ms_ and _scrm_ are used automatically.


Supported Summary Statistics
----------------------------
* sumstat_seg_sites: Segregating Sites
* sumstat_sfs: The Site Frequency Spectrum of a population.
* sumstat_jsfs: The Joint Site Frequency Spectrum of two populations.

Please refer to the statistic's help page for additional information (e.g.
`help(sumstat_jsfs)`).
