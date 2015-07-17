This fixes the invalid read found by valgrind:
http://www.stats.ox.ac.uk/pub/bdr/memtests/valgrind/coala/tests/testthat.Rout
Sorry for the inconvenience.

## Test environments
* Ubuntu 15.04 (local), R 3.2.1
* Ubuntu 12.04 (on Travis-CI), R 3.2.1
* Windows (on AppVeyor), R-devel
* win-builder (R-devel and R-release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTEs:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Paul Staab <develop@paulstaab.de>'

Days since last update: 3

License components with restrictions and base license permitting such:
  MIT + file LICENSE
File 'LICENSE':
  YEAR: 2015
  COPYRIGHT HOLDER: Paul Staab



* checking package dependencies ... NOTE
  No repository set, so cyclic dependency check skipped
  
  Seen only on win-builder.
  
