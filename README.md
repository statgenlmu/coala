coala
=====

[![Linux Build Status](https://travis-ci.org/statgenlmu/coala.png?branch=master)](https://travis-ci.org/statgenlmu/coala) 
[![Windows Build status](https://ci.appveyor.com/api/projects/status/uoduv0q64ddnqfva/branch/master?svg=true)](https://ci.appveyor.com/project/paulstaab/coala-02w83/branch/master)
[![Coverage Status](https://coveralls.io/repos/statgenlmu/coala/badge.svg?branch=master)](https://coveralls.io/r/statgenlmu/coala)
[![CRAN Status](http://www.r-pkg.org/badges/version/coala)](http://cran.r-project.org/web/packages/coala)

Coala is an R package that allows to simulate biological sequences according
to a given model of evolution.  It can call a number of efficient 
simulators based on
[coalescent theory](https://en.wikipedia.org/wiki/Coalescent_theory). 
All simulators can be combined with the program _seq-gen_ to simulate finite 
site mutation models. 
Coala also directly imports the simulation results into `R`, and allows to
calculate additional summary statistics from the results.


Installation
------------

The package can be installed from CRAN using

```R
install.packages("coala")
```

If you want to use the simulation programs `ms`, `msms` or `seqgen`, 
they need to be installed separately. This is described in the 
["Using External Simulators" vignette](https://cran.r-project.org/web/packages/coala/vignettes/coala-install.html).


Usage
-----
Coala comes with a
[vignette](https://cran.r-project.org/web/packages/coala/vignettes/coala-intro.html)
that explains the packages concepts and is a good place to start.
The [ABC vignette](https://cran.r-project.org/web/packages/coala/vignettes/coala-abc.html) gives an example
on how coala can be used to conduct the simulations for an Approximate Bayesian
Computation analysis.


Bugs
----
If you encounter problems when using _coala_, please 
[file a bug report](https://github.com/statgenlmu/coala/issues) or mail to
`develop (at) paulstaab.de`.


Supported Simulators
--------------------
The package supports the coalescent simulators _ms_, _scrm_ and _msms_.
All simulators can be combined with _seq-gen_ to simulate finite sites 
mutation models. The programs _ms_, _msms_ and _seq-gen_ must be installed 
manually on the system, while the R version of _scrm_ is used automatically.


Development
-----------
To follow or participate in the development of `coala`, please install the 
development version from GitHub using

```R
devtools::install_github('statgenlmu/coala')
```

on Linux and OS X. This requires that you have `devtools` and a compiler or 
Xcode installed. Bug reports and pull request on GitHub are highly appreciated.
The [extending coala vignette](https://cran.r-project.org/web/packages/coala/vignettes/coala-extend.html)
contains information on how to create new summary statistics and add simulators
to coala.
