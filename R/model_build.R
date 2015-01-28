#' Add a feature or parameter to a model
#' @export
#' @param e1 The Model to which the feature/parameter should be added
#' @param e2 The feature/parameter to add
#' @return The extended model
"+.CoalModel" <- function(e1, e2) {
  e2name <- deparse(substitute(e2)) # Passed throw for error messages
  addToModel(e2, e1, e2name)
}

addToModel <- function(x, model, x_name) UseMethod("addToModel")
addToModel.default <- function(x, model, x_name) {
  stop("Can not add `", x_name, "` to model")
}

addToModel.Parameter <- function(par, model, par_name) model

addToModel.Par_Range <- function(par, model, par_name) {
  if (par$get_name() %in% get_parameter_table(model))
    stop("There is already a parameter with name ", par.name)

  new_par <- data.frame(name=par$get_name(),
                        lower.range=par$get_range()[1],
                        upper.range=par$get_range()[2],
                        stringsAsFactors=F)

  model$parameters <- rbind(model$parameters, new_par)
  model
}

addToModel.Feature <- function(feat, model, feat_name) {
  model$features <- rbind(model$features, feat$get_table())
  for (parameter in feat$get_parameters()) {
    model <- model + parameter
  }
  if (feat$get_inter_locus_var()) {
    model <- addInterLocusVariation(model, feat$get_group())
  }
  model
}


# Deactived
# Seems you can't have +.DemographicModel and +.Feature and then add a
# feature to a model => create a common base class as in ggplot2
# "+.Feature" <- function(e1, e2) {
#   e2name <- deparse(substitute(e2)) # Passed throw for error messages
#   if (is.par_model(e2)) {
#     e1$add_parameter(e2)
#     return(e1)
#   } else if (is.feature(e2)) {
#     e1$add_feature(e2)
#     return(e1)
#   }
#   else stop("Can not add `", x_name, "` to feature")
# }
