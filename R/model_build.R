#' Add a feature or parameter to a model
#' @export
#' @param e1 The Model to which the feature/parameter should be added
#' @param e2 The feature/parameter to add
#' @return The extended model
#' @keywords internal
"+.coalmodelpart" <- function(e1, e2) {
  e2name <- deparse(substitute(e2)) # Passed through for error messages

  if (is.model(e1)) return(add_to_model(e2, e1, e2name))

  if (is_partial_model(e1)) {
    e1[[length(e1) + 1]] <- e2
    return(e1)
  }

  partial_model <- list(e1, e2)
  class(partial_model) <- c("partial_model", "coalmodelpart")
  partial_model
}

is_partial_model <- function(x) inherits(x, "partial_model")

#' @export
print.partial_model <- function(x, ...) {
  cat("Partial coal_model with ", length(x), " components:\n")
  for (part in x) {
    cat("* ")
    print(part)
  }
}

add_to_model <- function(x, model, x_name) UseMethod("add_to_model")

#' @export
add_to_model.default <- function(x, model, x_name) {
  stop("Can not add `", x_name, "` to model")
}

#' @export
add_to_model.parameter <- function(x, model, x_name) model


#' @export
add_to_model.named_par <- function(x, model, x_name) {
  if (x$get_name() %in% get_par_names(model))
    stop("There is already a parameter with name ", x_name)

  model$parameter[[length(model$parameter) + 1]] <- x

  model$id <- get_id()
  model
}


#' @export
add_to_model.variation_par <- function(x, model, x_name) {
  model <- add_variation(model)
  for (par in x$get_base_par()) model <- model + par
  model$id <- get_id()
  model
}


#' @export
add_to_model.feature <- function(x, model, x_name) {
  # Check that the population in the feature exists
  pop <- x$get_population()
  if (!is.null(pop)) {
    if (any(!pop %in% get_populations(model)) && pop[1] != "all" ) {
      stop("Invalid population in ", x_name, call. = FALSE)
    }
  }

  # Check that the locus group is set correctly
  locus_group <- x$get_locus_group()
  if (is.null(locus_group)) {
    stop(x_name, ": locus group has an invalid value", call. = FALSE)
  }

  # Execute the features checks
  x$check(model)

  # Add the parameters in the feature
  for (para in x$get_parameters()) model <- model + para

  # Add the feature itself
  model$features[[length(model$features) + 1]] <- x

  model$id <- get_id()
  model
}


#' @export
add_to_model.locus <- function(x, model, x_name) {
  model$loci[[length(model$loci) + 1]] <- x
  model$id <- get_id()
  model
}


#' @export
add_to_model.partial_model <- function(x, model, x_name) {
  for (part in x) {
    model <- model + part
  }
  model
}
