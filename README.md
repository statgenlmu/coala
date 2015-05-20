coala
=====

[![Linux Build Status](https://travis-ci.org/statgenlmu/coala.png?branch=master)](https://travis-ci.org/statgenlmu/coala) 
[![Windows Build status](https://ci.appveyor.com/api/projects/status/uoduv0q64ddnqfva/branch/master?svg=true)](https://ci.appveyor.com/project/paulstaab/coala-02w83/branch/master)
[![Coverage Status](https://coveralls.io/repos/statgenlmu/coala/badge.svg?branch=master)](https://coveralls.io/r/statgenlmu/coala)

This is a framework for calling coalescent simulators from within R. It allows to 
specify a model using a ggplot2-like syntax and supports multiple simulators. It can
also incorporate the program _seq-gen_ to simulate finite site mutation models.

This package is currently in early development. An initial CRAN release will occur when the
milestone `0.1` is reached.

The development version can be installed via

```R
devtools::install_github('statgenlmu/coala')
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

