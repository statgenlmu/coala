#' @export
print.Coalmodel <- function(x, ...) {
  if (length(get_features(x)) == 0) cat("The model has no features\n")
  else {
    cat("Features:\n")
    for (feat in get_features(x)) {
      cat("* ")
      feat$print()
    }
  }

  if (length(get_parameter(x)) == 0) cat("The model has no parameters\n")
  else {
    cat("\nParameters:\n")
    for (par in get_parameter(x)) {
      cat("* ")
      par$print()
    }
  }

  if (length(get_loci(x)) == 0) cat("The model has no loci\n")
  else {
    cat("\nLoci:\n")
    for (locus in get_loci(x)) {
      cat("* ")
      locus$print()
    }
    cat("\n")
  }

  cat("\nSimulator:", select_simprog(x)$get_name(), "\n")
  cmd <- "Can not simulate model"
  tryCatch(cmd <- get_cmd(x), error = function(e) NULL)
  cat("Command:", paste(cmd, collapse="\n"), "\n")
}
