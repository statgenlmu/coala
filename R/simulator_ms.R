#' Generate command line arguments for features
#'
#' These functions are exported only for technical reasons
#' (because they are S3 methods) and are not indended for
#' users.
#'
#' @param feature The feature for which the argument is generated
#' @param model The complete model for which the argument is generated
conv_to_ms_arg <- function(feature, model) UseMethod("conv_to_ms_arg")

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_ms_arg.default <- function(feature, model) {
  stop("Unknown feature", call. = FALSE)
}


ms_create_cmd_tempalte <- function(model) {
  cmd <- read_cache(model, 'ms_cmd')
  if (is.null(cmd)) {
    cmd <- paste(vapply(model$features, conv_to_ms_arg,
                        FUN.VALUE = character(1), model),
                 collapse = "")
    cmd <- paste0("c('", cmd, "')")
    cache(model, 'ms_cmd', cmd)
  }

  cmd
}


#' @importFrom phyclust ms
#' @importFrom R6 R6Class
#' @include simulator_class.R
Simulator_ms <- R6Class('Simulator_ms', inherit = Simulator,
  private = list(
    name = 'ms',
    priority = 100
  ),
  public = list(
    simulate = function(model, parameters=numeric()) {
      stopifnot(length(parameters) == 0 | all(is.numeric(parameters)))

      template <- ms_create_cmd_tempalte(model)
      sample_size <- sum(get_sample_size(model, for_sim = TRUE))

      # Do the actuall simulation
      files <- lapply(1:get_locus_group_number(model) , function(i) {
        sim_cmds <- fill_cmd_template(template, model, parameters, i)
        #print(sim_cmds)


        vapply(1:nrow(sim_cmds), function(i) {
          file <- tempfile("ms")
          ms(sample_size,
             format(sim_cmds[i, "locus_number"], scientific = FALSE),
             sim_cmds[i, "command"], file)
          file
        }, character(1))
      })

      # Parse the output and calculate summary statistics
      if (requires_segsites(model)) {
        seg_sites <- parse_ms_output(files,
                                     get_sample_size(model, for_sim = TRUE),
                                     get_locus_number(model))

        if (has_trios(model)) {
          seg_sites <- conv_for_trios(seg_sites, model)
        }
      } else {
        seg_sites <- NULL
      }

      sum_stats <- calc_sumstats(seg_sites, files, model, parameters)

      # Clean Up
      unlink(unlist(files))
      sum_stats
    },
    get_cmd = function(model) {
      template <- ms_create_cmd_tempalte(model)
      cmd <- fill_cmd_template(template, model, NULL, 1, eval_pars = FALSE)
      paste("ms",
            sum(get_sample_size(model, TRUE)),
            cmd[1, "locus_number"],
            cmd[1, "command"])
    }
  )
)

register_simulator(Simulator_ms)
