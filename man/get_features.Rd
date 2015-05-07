% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/feature_sample.R, R/model_getters.R, R/sumstat.R
\name{get_sample_size}
\alias{get_features}
\alias{get_locus_length}
\alias{get_locus_length_matrix}
\alias{get_locus_number}
\alias{get_outgroup}
\alias{get_outgroup_size}
\alias{get_parameter}
\alias{get_parameter_table}
\alias{get_population_indiviuals}
\alias{get_populations}
\alias{get_sample_size}
\alias{get_summary_statistics}
\title{Getters for coalescent models}
\usage{
get_sample_size(model, for_sim = FALSE)

get_features(model)

get_parameter_table(model)

get_parameter(model)

get_locus_length(model, locus = NULL, group = NULL, total = TRUE)

get_populations(model)

get_locus_length_matrix(model)

get_locus_number(model, group = NA)

get_outgroup(model)

get_outgroup_size(model, for_sim = FALSE)

get_population_indiviuals(model, pop, zero_indexed = FALSE)

get_summary_statistics(model)
}
\arguments{
\item{model}{The coalescent model from which aspects are returned}

\item{for_sim}{If true, the sample size used internally for the simulation
will be reported rather than the number of actuall samples. The numbers
can be unequal for the simulation of unphased data.}

\item{locus}{The number of the locus.}

\item{total}{If \code{FALSE}, the length of loci in a trio will be reported
  individually. If \code{TRUE} the sum of the loci's length will be reported.
  This does not affect non-trio loci.}

\item{pop}{The population for which aspects are returned}
}
\description{
Getters for coalescent models
}
\section{Functions}{
\itemize{
\item \code{get_sample_size}: Returns a vector of samples sizes per
population.

\item \code{get_parameter_table}: Returns the ranged parameters of a model as a
data.frame

\item \code{get_parameter}: Returns the ranged parameters of a model

\item \code{get_locus_length}: Returns the length of the loci in a locus group

\item \code{get_populations}: Returns a vector of populations in the model

\item \code{get_locus_length_matrix}: Returns a matrix with detailed length
information about the loci in the model.

\item \code{get_locus_number}: Returns the number of loci in a locus group

\item \code{get_outgroup}: Returns the population that is marked as outgroup

\item \code{get_outgroup_size}: Returns the number of samples in the outgroup

\item \code{get_population_indiviuals}: Returns the index of the individuals of one
population

\item \code{get_summary_statistics}: Returns the summary statistics in the model
}}
\author{
Paul Staab
}
