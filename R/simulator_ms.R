# Translating the model into simulation commands
conv_to_ms_arg <- function(feature, model) UseMethod("conv_to_ms_arg")
conv_to_ms_arg.default <- function(feature, model) {
  stop("Unknown feature when generating ms command")
}


# This function generates an string that contains an R command for generating
# an ms call to the current model.
ms_generate_opts_cmd <- function(model) {
  cmd <- paste(vapply(model$features, conv_to_ms_arg,
                      FUN.VALUE = character(1), model),
               collapse = "")
  paste0("c('", cmd, "')")
}


ms_generate_opts <- function(model, parameters, group, eval_pars = TRUE) {
  if (eval_pars) locus_length <- get_locus_length(model, group = group)
  else locus_length <- "locus_length"

  ms_tmp <- create_par_env(model, parameters, locus = group,
                           locus_length = locus_length,
                           for_cmd = !eval_pars)

  cmd <- read_cache(model, 'ms_cmd')
  if (is.null(cmd)) {
    cmd <- ms_generate_opts_cmd(model)
    cache(model, 'ms_cmd', cmd)
  }

  eval(parse(text = cmd), envir = ms_tmp)
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

      # Do the actuall simulation
      files <- lapply(1:get_locus_group_number(model) , function(i) {
        opts <- ms_generate_opts(model, parameters, i)
        #print(opts)
        file <- tempfile('csr_ms')

        ms(sum(get_sample_size(model, for_sim = TRUE)),
           format(get_locus_number(model, group = i), scientific = FALSE),
           unlist(strsplit(opts, " ")), file)

        if (file.info(file)$size == 0) stop("ms simulation output is empty")
        file
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
      unlink(files)
      sum_stats
    },
    get_cmd = function(model) {
      cmd <- ms_generate_opts(model, NULL, 1, FALSE)
      paste("ms",
            sum(get_sample_size(model, TRUE)),
            get_locus_number(model),
            paste(cmd, collapse = ' '))
    }
  )
)

register_simulator(Simulator_ms)
