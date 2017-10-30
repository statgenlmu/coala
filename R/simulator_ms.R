#' Generate command line arguments for features
#'
#' These functions are exported only for technical reasons
#' (because they are S3 methods) and are not intended for
#' users.
#'
#' @param feature The feature for which the argument is generated
#' @param model The complete model for which the argument is generated
#' @keywords internal
conv_to_ms_arg <- function(feature, model) UseMethod("conv_to_ms_arg")

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_ms_arg.default <- function(feature, model) {
  stop("Unknown feature", call. = FALSE)
}


#' @importFrom R6 R6Class
#' @include simulator_class.R
ms_class <- R6Class("ms", inherit = simulator_class,
  private = list(
    name = "ms",
    priority = 100,
    binary = NULL
  ),
  public = list(
    initialize = function(priority = 300) {
      if (!requireNamespace("phyclust", quietly = TRUE)) {
        stop("Please install package 'phyclust' to use ms", call. = FALSE)
      }
      assert_that(is.numeric(priority) && length(priority) == 1)
      private$priority <- priority
    },
    create_cmd_template = function(model) {
      cmd <- read_cache(model, "ms_cmd")
      if (is.null(cmd)) {
        cmd <- paste(vapply(model$features, conv_to_ms_arg,
                            FUN.VALUE = character(1), model),
                     collapse = "")
        cmd <- paste0("c('", cmd, "')")
        cache(model, "ms_cmd", cmd)
      }

      cmd
    },
    create_task = function(model, pars, locus_number,
                           locus_id = 1,
                           eval_pars = TRUE) {
      tmplt <- self$create_cmd_template(model)
      cmd <- fill_cmd_template(tmplt, model, pars, locus_id, eval_pars)
      create_sim_task(self, locus_number,
                      cmd = cmd,
                      sample_size = sum(get_sample_size(model, for_sim = TRUE)))
    },
    call = function(sample_size, n_loci, command) {
      output <- tempfile("ms")
      phyclust::ms(sample_size, n_loci, command, temp.file = output)
      if (!file.exists(output)) stop("ms simulation failed")
      list(file = output, cmd = paste("ms", sample_size, n_loci, command))
    },
    simulate = function(model, sim_task) {
      # Call ms
      result <- self$call(sim_task$get_arg("sample_size"),
                          sim_task$locus_number,
                          sim_task$get_arg("cmd"))

      # Parse the output
      if (requires_segsites(model) || requires_trees(model)) {
        output <- parse_ms_output(list(result$file),
                                  get_sample_size(model, for_sim = TRUE),
                                  sim_task$locus_number)
      } else {
        output <- list(seg_sites = NULL, trees = NULL)
      }

      # Add the file if needed
      if (requires_files(model)) output$files <- result$file
      else unlink(result$file)

      # Add the simulation cmd
      output$cmd <- result$cmd
      output$simulator <- self

      output
    },
    get_cmd = function(model) {
      task <- self$create_task(model, NULL, get_locus_number(model), 1, FALSE)
      paste(self$get_name(),
            task$get_arg("sample_size"),
            task$locus_number,
            task$get_arg("cmd"))
    },
    get_info = function() c(name = "ms",
                            version = paste0("phyclust_",
                                             packageVersion("phyclust")))
  )
)

has_ms <- function() !is.null(simulators[["ms"]])


#' Simulator: ms
#'
#' This function adds the simulator 'ms' to the list of available simulators.
#' In order to use 'ms', you need to install the CRAN package \pkg{phyclust}.
#' By default, 'scrm' will be preferred over 'ms'. Raise the priority of 'ms'
#' to change this behavior.
#'
#' @references
#' Richard R. Hudson.
#' Generating samples under a Wright-Fisher neutral model of genetic variation.
#' Bioinformatics (2002) 18 (2): 337-338
#' doi:10.1093/bioinformatics/18.2.337
#'
#' @name simulator_ms
#' @param priority The priority for this simulator. If multiple simulators
#'   can simulate a model, the one with the highest priority will be used.
#' @export
#' @examples
#' # To prefer ms to scrm:
#' \dontrun{activate_ms(priority = 500)}
#' @family simulators
activate_ms <- function(priority = 300) {
  register_simulator(ms_class$new(priority))
  reset_cache()
  invisible(NULL)
}
