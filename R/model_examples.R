# Example models used in unittests

model_theta_tau <- function() {
  dm.createDemographicModel(c(10, 15), 10) +
    feat_pop_merge(par_range('tau', 0.01, 5), 2, 1) +
    feat_mutation(par_range('theta', 1, 10))
}
