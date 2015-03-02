#' Add a feature or parameter to a model
#' @export
#' @param e1 The Model to which the feature/parameter should be added
#' @param e2 The feature/parameter to add
#' @return The extended model
"+.Coalmodel" <- function(e1, e2) {
  e2name <- deparse(substitute(e2)) # Passed throw for error messages
  add_to_model(e2, e1, e2name)
}


add_to_model <- function(x, model, x_name) UseMethod("add_to_model")
add_to_model.default <- function(x, model, x_name) {
  stop("Can not add `", x_name, "` to model")
}


add_to_model.Coalmodel <- function(add, model, par_name) {
  for (feat in get_features(add)) model <- model + feat
  for (par in get_parameter(add)) model <- model + par
  for (stat in add$sum_stats) model <- model + stat

  model$loci <- rbind(model$loci, add$loci)
  model
}


add_to_model.Parameter <- function(par, model, par_name) model


add_to_model.Par_Range <- function(par, model, par_name) {
  if (par$get_name() %in% get_par_names(model))
    stop("There is already a parameter with name ", par_name)

  model$parameter[[length(model$parameter) + 1]] <- par

  model$id <- get_id()
  model
}


add_to_model.Feature <- function(feat, model, feat_name) {
  for (para in feat$get_parameters()) model <- model + para
  feat$reset_parameters()

  model$features[[length(model$features) + 1]] <- feat

  if (feat$get_inter_locus_var()) {
    model <- add_inter_locus_var(model, feat$get_group())
  }

  model$id <- get_id()
  model
}


add_to_model.Locus <- function(locus, model, locus_name) {
  if (length(locus$get_name()) == 1) {
    locus_names <- c('', locus$get_name(), '')
  } else if (length(locus$get_name()) == 3) {
    locus_names <- locus$get_name()
  } else stop("Failed to get locus names from ", locus_name)

  if (length(locus$get_length()) == 1) {
    locus_length <- c(0, 0, locus$get_length(), 0, 0)
  } else if (length(locus$get_length()) == 5) {
    locus_length <- locus$get_length()
  } else stop("Failed to get locus length from ", locus_name)

  model$loci <- rbind(model$loci, data.frame(group = locus$get_group(),
                                             number = locus$get_number(),
                                             name = locus_names[2],
                                             name_l = locus_names[1],
                                             name_r = locus_names[3],
                                             length_l = locus_length[1],
                                             length_il = locus_length[2],
                                             length_m = locus_length[3],
                                             length_ir = locus_length[4],
                                             length_r = locus_length[5],
                                             stringsAsFactors = FALSE))

  model$id <- get_id()
  model
}
