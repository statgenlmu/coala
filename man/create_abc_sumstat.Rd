% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interface_abc.R
\name{create_abc_sumstat}
\alias{create_abc_sumstat}
\title{Convert Simulation Results to abc's Summary Statistic Format}
\usage{
create_abc_sumstat(sim_results, model)
}
\arguments{
\item{sim_results}{The simulation results as returned from
\code{\link[=simulate.coalmodel]{simulate}}.}

\item{model}{The model used for the simulations.}
}
\value{
A data.frame that can be used as \code{sumstat}
argument of \code{\link[abc]{abc}}.
}
\description{
This function creates an object compatible with the \code{sumstat}
argument of the \code{\link[abc]{abc}} function from coala's simulation
results. It converts all summary statistics that are in the simulation
results and expects that each of them is a numeric vector.
Use transformation functions to convert none vector-valued  statistics
(e.g. \code{\link{sumstat_jsfs}}, \code{\link{sumstat_omega}} or
\code{\link{sumstat_trees}}) into a vector.
}
\examples{
# Using the SFS:
model <- coal_model(10, 1) +
  feat_mutation(par_prior("theta", rnorm(1, 5, .5))) +
  sumstat_sfs()
sim_results <- simulate(model, nsim = 2)
abc_sumstat <- create_abc_sumstat(sim_results, model)
print(abc_sumstat)

# Using the JSFS and converting it into a vector:
model <- coal_model(c(10, 10), 1) +
  feat_mutation(par_prior("theta", rnorm(1, 5, .5))) +
  feat_migration(par_prior("m", rnorm(1, .5, .1)), symmetri = TRUE) +
  sumstat_jsfs(transformation = function(jsfs) {
    c(sum(jsfs[1, ]), sum(jsfs[, 1]), sum(jsfs[-1, -1]))
  })
sim_results <- simulate(model, nsim = 2)
abc_sumstat <- create_abc_sumstat(sim_results, model)
print(abc_sumstat)
}
\seealso{
For generating abc's parameter format:
\code{\link{create_abc_param}}
}
