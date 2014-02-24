.PHONY: howtos install test test-setup integration-test travis-test check clean

VERSION=$(shell grep Version DESCRIPTION | awk '{print $$2}')
PACKAGE=EvolutionaryModel_$(VERSION).tar.gz
R_CHECK_ARGS?="--as-cran"
R_BUILD_ARGS?=

R_SOURCES=$(wildcard R/*.R) 
CPP_SOURCES=$(wildcard src/*.cc)
VIGNETTES=$(wildcard vignettes/*.pdf)
TESTS=$(wildcard inst/unitTests/*.R) $(wildcard tests/*.R)

default: $(PACKAGE)

release: clean test-setup howtos $(PACKAGE) check  
travis-test: $(PACKAGE) test-setup 

test: install
	cd tests; Rscript --vanilla test-all.r

check: $(PACKAGE)
	# Runs an R CMD check
	R CMD check $(R_CHECK_ARGS) $(PACKAGE)

package: $(PACKAGE) 

install: 
	R CMD INSTALL .

$(PACKAGE): $(R_SOURCES) $(CPP_SOURCES) $(TESTS) $(VIGNETTES) DESCRIPTION man
	R CMD build $(R_BUILD_ARGS) .

README: README.md
	grep -v "\`\`\`" README.md | grep -v "Build Status" > README

man: $(R_SOURCES) DESCRIPTION
	- rm -r man 2> /dev/null
	Rscript -e 'library(roxygen2); roxygenise(".")'

clean:
	- rm -rv jaatha.Rcheck
	- cd src/; rm *.so *.o *.rds ms/*.o 2> /dev/null
	- rm -r man 2> /dev/null
	- cd howtos; make clean
	- rm -rv inst/unitTests/test_setup.Rda
