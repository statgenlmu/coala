% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_growth.R
\name{feat_growth}
\alias{feat_growth}
\title{Feature: Exponential population size growth/decline}
\usage{
feat_growth(rate, population = "all", time = "0", locus_group = "all")
}
\arguments{
\item{rate}{The growth rate. Can be a numeric or a \code{\link{parameter}}.
See \code{Details} for an explanation how the rate affects the
population size.}

\item{population}{The population which growths/declines. Can be
"all" for all populations, or the number of one population.}

\item{time}{The time at which the growth rate is changed. Can also be
a \code{\link{parameter}}.}

\item{locus_group}{The loci for which this features is used. Can either be
\code{"all"} (default), in which case the feature is used for simulating
all loci, or a numeric vector. In the latter case, the feature is only
used for the loci added in \code{locus_} commands  with the corresponding
index starting from 1 in order in which the commands where added to the
model. For example, if a model has
\code{locus_single(10) + locus_averaged(10, 11) + locus_single(12)} and
this argument is \code{c(2, 3)}, than the feature is used for all but
the first locus (that is locus 2 - 12).}
}
\value{
The feature, which can be added to a model created with
\code{\link{coal_model}} using \code{+}.
}
\description{
This feature changes the growth factor of a population at given
point in time. This factor applies to the time interval further
into the past from this point.
}
\details{
The population size changes by a factor \eqn{exp(-\alpha*t)}, where
\eqn{\alpha} is the growth parameter and \eqn{t} is the time since
the growth has started. For positive alpha, the population will decline
backwards in time or grow forwards in time. For a negative value of
\eqn{\alpha} it will decline (forward in time).
}
\examples{
# Simulate a haploid population that has been expanding for
# the last 2*Ne generations
model <- coal_model(10, 1) +
  feat_growth(5, time = 0) +
  feat_growth(0, time = 1) +
  feat_mutation(10) +
  sumstat_sfs()
simulate(model)
}
\seealso{
For instantaneous population size
changes: \code{\link{feat_size_change}}

For creating a model: \code{\link{coal_model}}

Other features: 
\code{\link{feat_ignore_singletons}()},
\code{\link{feat_migration}()},
\code{\link{feat_mutation}()},
\code{\link{feat_outgroup}()},
\code{\link{feat_pop_merge}()},
\code{\link{feat_recombination}()},
\code{\link{feat_selection}()},
\code{\link{feat_size_change}()},
\code{\link{feat_unphased}()}
}
\concept{features}
