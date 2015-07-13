#' Creates a coalescent model
#'
#' This creates a basic coalescent model to which more features, loci,
#' parameters and summary statistics can be added. Data under the model
#' can be simulated using the \code{simulate} function.
#'
#' @param sample_size Defines the number of populations and the number of
#'   samples taken from each population. Given as an integer vector which
#'   length defines the number of population and each entry gives the number
#'   of haploids sampled from the corresponding population.
#' @param loci_number You can optionally add a number of loci with equal length
#'   to the model. This gives to number of loci to add.
#' @param loci_length This gives the length of the loci to add.
#' @return The coalescent model
#' @export
#' @importFrom assertthat assert_that
#' @examples
#' model <- coal_model(5, 10, 100) + feat_mutation(1) + sumstat_sfs()
#' simulate(model)
#'
#' model <- coal_model(c(2, 2), 1, 10) +
#'   feat_recombination(1) +
#'   feat_migration(0.5, symmetric = TRUE) +
#'   sumstat_trees()
#' simulate(model)
coal_model <- function(sample_size, loci_number=0, loci_length=1000) {
  model <- list()
  class(model) <- c("coalmodel")

  model$features <- list()
  model$loci <- list()
  model$parameter <- list()
  model$sum_stats <- create_sumstat_container()
  model$has_variation <- FALSE

  model$scaling_factor <- 1
  model$id <- get_id()

  # Add sample sizes
  assert_that(is.numeric(sample_size))
  assert_that(sum(sample_size) > 0)
  model <- model + feat_sample(sample_size)

  # Add locus
  if (loci_number > 0) {
    model <- model + locus_averaged(loci_number, loci_length)
  }

  model
}


is.model <- function(model) {
  "coalmodel" %in% class(model)
}


# Selects a program for simulation that is capable of all current features
select_simprog <- function(model) {
  name <- read_cache(model, "simprog")

  if (is.null(name)) {
    priority <- -Inf

    for (simprog_name in ls(simulators)) {
      simprog <- get_simulator(simprog_name)
      valid <- try(simprog$get_cmd(model), silent = TRUE)
      if (all(class(valid) != "try-error")) {
        if (simprog$get_priority() > priority) {
          name <- simprog
          priority <- simprog$get_priority()
        }
      }
    }

    if (is.null(name)) warning("No suitable simulation software found!")
    cache(model, "simprog", name)
  }

  name
}


add_variation <- function(model) {
  model$has_variation <- TRUE
  model
}


has_variation <- function(model) model$has_variation


has_trios <- function(model) {
  sum(get_locus_length_matrix(model)[ , c(1:2, 4:5)]) > 0
}


get_snp_positions <- function(seg_sites, model, relative=TRUE) {
  lapply(1:length(seg_sites), function(locus) {
    pos <- attr(seg_sites[[locus]], "position")
    locus_length <- get_locus_length(model, locus, total = FALSE)

    # Nothing changes without trios
    if (length(locus_length) == 1) {
      if (relative) return(pos)
      else return(pos * locus_length)
    }

    # Convert if we have trios
    trio_locus <- attr(seg_sites[[locus]], "locus")
    if (is.null(trio_locus)) trio_locus <- 0
    pos[trio_locus == -1] <- pos[trio_locus == -1] * locus_length[1]
    pos[trio_locus == 0] <- pos[trio_locus == 0] * locus_length[3] +
      sum(locus_length[1:2])
    pos[trio_locus == 1] <- pos[trio_locus == 1] * locus_length[5] +
      sum(locus_length[1:4])
    if (relative) pos <- pos / sum(locus_length)
    pos
  })
}
