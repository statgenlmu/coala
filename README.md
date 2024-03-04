coala
=====

Coala is an R package for simulating biological sequences according
to a given model of evolution.  It can call a number of efficient 
simulators based on
[coalescent theory](https://en.wikipedia.org/wiki/Coalescent_theory). 
All simulators can be combined with the program _seq-gen_ to simulate finite 
site mutation models. 
Coala also directly imports the simulation results into `R`, and can
calculate various summary statistics from the results.


Installation
------------

The package can be installed from CRAN using

```R
install.packages("coala")
```

If you want to use the simulation programs `ms`, `msms` or `seqgen`, 
they need to be installed separately. This is described in the 
["Using External Simulators" vignette](https://github.com/statgenlmu/coala/blob/master/vignettes/coala-install.Rmd) and
in [the wiki](https://github.com/statgenlmu/coala/wiki/Installation).


Usage & Help
------------
Coala comes with a
[vignette](https://github.com/statgenlmu/coala/blob/master/vignettes/coala-intro.Rmd)
that explains the packages concepts and is a good place to start. It also has a 
[vignette containing a few example applications](https://github.com/statgenlmu/coala/blob/master/vignettes/coala-examples.Rmd).

Detailed information about coala's functions is provided via R's help system. 
Call `help(_function_)` in R to view them. They usually also contain examples and further links.

The [ABC vignette](https://github.com/statgenlmu/coala/blob/master/vignettes/coala-abc.Rmd) 
gives an example on how coala can be used to conduct the simulations for [Approximate Bayesian
Computation](https://en.wikipedia.org/wiki/Approximate_Bayesian_computation).

Also take a look at the [project wiki](https://github.com/statgenlmu/coala/wiki) for additional
resources.


Example
-------
In the following example, we create a simple panmictic model, simulate it and 
calculate the site frequency spectrum (SFS) of the simulation results:

```R
model <- coal_model(sample_size = 10, loci_number = 2) +
  feat_mutation(5) +
  sumstat_sfs()
result <- simulate(model)
result$sfs
# [1] 15 12  1  4  0  1  0  2  0
```

More examples can be found in the 
[examples vignette](https://github.com/statgenlmu/coala/blob/master/vignettes/coala-examples.Rmd).


Problems
--------
If you encounter problems when using _coala_, please 
[file a bug report](https://github.com/statgenlmu/coala/issues) or mail to
`coala-pkg (at) googlegroups.com`.


Supported Simulators
--------------------
The package supports the coalescent simulators _ms_, _scrm_ and _msms_.
All simulators can be combined with _seq-gen_ to simulate finite sites 
mutation models. The programs _msms_ and _seq-gen_ must be [installed 
manually](https://github.com/statgenlmu/coala/wiki/Installation#installing-additional-simulators). 
The R version of _scrm_ should be installed automatically,
and the R version _ms_ if the package `phyclust` is installed.


Development
-----------
To follow or participate in the development of `coala`, please install the 
development version from GitHub using

```R
devtools::install_github('statgenlmu/coala')
```

on Linux and OS X. This requires that you have `devtools` and a compiler or 
Xcode installed. Bug reports and pull request on GitHub are highly appreciated.
The [extending coala vignette](https://github.com/statgenlmu/coala/blob/master/vignettes/coala-extend.Rmd)
contains information on how to create new summary statistics and add simulators
to coala. The [wiki](https://github.com/statgenlmu/coala/wiki) also contains a few
resources for developers.
