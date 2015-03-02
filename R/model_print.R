#' @export
print.Coalmodel <- function(x, ...) {
  if (length(get_features(x)) == 0) cat("The model has no features\n")
  else {
    cat("Features:\n")
    for (feat in get_features(x)) print(feat)
  }

  if (length(get_parameter(x)) == 0) cat("The model has no parameters\n")
  else {
    cat("Parameters:\n")
    for (par in get_parameter(x)) print(par)
  }
}
