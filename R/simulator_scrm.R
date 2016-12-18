conv_to_scrm_arg <- function(feature, model) UseMethod("conv_to_scrm_arg")

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_scrm_arg.default <- function(feature, model) {
  stop("Unknown feature when generating scrm command")
}


#' @importFrom scrm scrm
#' @include simulator_class.R
#' @include simulator_ms.R
scrm_class <- R6Class('Scrm', inherit = ms_class, #nolint
  private = list(
    name = "scrm",
    version = packageDescription("scrm", fields = "Version")
  ),
  public = list(
    create_cmd_template = function(model) {
      cmd <- read_cache(model, "scrm_cmd")
      if (is.null(cmd)) {
        cmd <- paste(vapply(model$features, conv_to_scrm_arg,
                            FUN.VALUE = character(1), model),
                     collapse = "")
        cmd <- paste0("c('", cmd, "')")
        cache(model, "scrm_cmd", cmd)
      }
      cmd
    },
    initialize = function(priority = 400) {
      assert_that(is.numeric(priority) && length(priority) == 1)
      private$priority <- priority
    },
    simulate = function(model, sim_task) {
      sample_size <- sum(get_sample_size(model, for_sim = TRUE))

      if (requires_files(model)) file <- tempfile("scrm")
      else file <- ""

      cmd <- paste(sample_size, sim_task$locus_number, sim_task$get_arg("cmd"))
      result <- scrm(cmd, file)

      if (requires_segsites(model)) {
        result$seg_sites <- lapply(result$seg_sites, function(x) {
          create_segsites(x, as.numeric(colnames(x)), check = FALSE)
        })
      }

      if (requires_trees(model) && packageVersion("scrm") < "1.7.2") {
        result$trees <- lapply(result$trees, function(x) {
          strsplit(x, split = "\n", fixed = TRUE)[[1]]
        })
      }

      result$cmd <- paste("scrm", cmd)
      result$simulator <- self

      if (file != "") result$files <- file

      result
    },
    get_info = function() c(name = "scrm", version = private$version)
  )
)

#' Simulator: scrm
#'
#' This function adds the simulator 'scrm' to the list of available simulators.
#' It is provided via the CRAN package \pkg{scrm} and should be always installed
#' alongside with \pkg{coala}. It should be activated automatically, and this
#' function is only needed to change it \code{priority}.
#'
#' @references
#' Paul R. Staab, Sha Zhu, Dirk Metzler and Gerton Lunter (2015).
#' "scrm: efficiently simulating long sequences using the approximated
#' coalescent with recombination."
#' Bioinformatics, 31(10), pp. 1680-1682.
#' http://dx.doi.org/10.1093/bioinformatics/btu861
#'
#' @name simulator_scrm
#' @inheritParams simulator_ms
#' @family simulators
#' @export
#' @examples
#' # Change scrm's priority
#' model <- coal_model(10, 1) + feat_mutation(5)
#' model # scrm is used by default
#' activate_scrm(250)
#' model # Now ms is used instead (if installed)
#' activate_scrm(550)
#' model # Now scrm is used again
activate_scrm <- function(priority = 400) {
  register_simulator(scrm_class$new(priority))
  reset_cache()
  invisible(NULL)
}
