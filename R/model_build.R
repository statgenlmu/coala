#' Add a feature or parameter to a model
#' @export
#' @param e1 The Model to which the feature/parameter should be added
#' @param e2 The feature/parameter to add
#' @return The extended model
"+.DemographicModel" <- function(e1, e2) {
  e2name <- deparse(substitute(e2)) # Passed throw for error messages
  addToModel(e2, e1, e2name)
}

addToModel <- function(x, model, x_name) UseMethod("addToModel")
addToModel.default <- function(x, model, x_name) {
  stop("Can not add `", x_name, "` to model")
}

addToModel.Parameter <- function(par, model, par_name) model

addToModel.Par_Range <- function(par, model, par_name) {
  range <- par$get_range()
  dm.addParameter(model, par$get_name(), range[1], range[2])
}

addToModel.Feature <- function(feat, model, feat_name) {
  dm.addFeature(model, feat)
}