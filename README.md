coalescentsimulatr
=================

I developed an interface to call coalescent simulation programs from within R for 
my CRAN package [jaatha](https://github.com/paulstaab/jaatha). To make this 
functionality available to an broader audience (and to clean up jaatha), I aim
to separate this interface into a new package. Along the way, I plan to polish
the somewhat clumsy syntax and to switch from S4 classes to Reference Classes or
to the new R6 system.

Coalescent Simulation
---------------------
Coalescent/backwards simulation allow to rapidly simulated genetic data according 
to a given model of evolution. The key idea is to trace the ancestry of the
sampled individuals/chromosomes backwards in time using the mathematical framework of 
coalescence theory, rather than simulation the evolution of the complete population
forwards in time. This is several orders of magnitude faster, at least as long 
as the assumptions of coalescence theory hold (population is large compared to 
the sample, random mating, Wright-Fisher-dynamics, ...).

Supported Simulators
--------------------
There is a large number of coalescent simulator available. The code in jaatha 
currently supports simulations using

* ms
* ms + seq-gen
* msms

where seq-gen and msms must be installed manually.
