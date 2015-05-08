# Translating the model into simulation commands
conv_to_scrm_arg <- function(feature, model) UseMethod("conv_to_scrm_arg")
conv_to_scrm_arg.default <- function(feature, model) {
  stop("Unknown feature when generating scrm command")
}


scrm_generate_opts_cmd <- function(model) {
  cmd <- read_cache(model, "scrm_cmd")
  if (is.null(cmd)) {
    cmd <- paste(vapply(model$features, conv_to_scrm_arg,
                        FUN.VALUE = character(1), model),
                 collapse = "")
    cmd <- paste0("c(",
                  sum(get_sample_size(model, for_sim = TRUE)), ", ",
                  "format(locus_number, scientific = FALSE), '",
                  cmd, "')")
    cache(model, "scrm_cmd", cmd)
  }
  cmd
}


scrm_generate_opts <- function(model, parameters, group, eval_pars = TRUE) {
  locus_length <- get_locus_length(model, group = group)
  locus_number <- get_locus_number(model, group = group)

  par_env <- create_par_env(model, parameters, locus = group,
                            locus_length = locus_length,
                            locus_number = locus_number,
                            for_cmd = !eval_pars)

  cmd <- scrm_generate_opts_cmd(model)

  eval(parse(text = cmd), envir = par_env)
}


#' @importFrom scrm scrm
#' @include simulator_class.R
#' @include simulator_ms.R
SimulatorScrm <- R6Class('SimulatorScrm', inherit = Simulator, #nolint
  private = list(
    name = 'scrm',
    priority = 90
  ),
  public = list(
    simulate = function(model, parameters) {
      if (length(get_locus_group_number(model)) > 1)
        stop("scrm can only simulate one group of loci at the moment.")

      args <- paste(scrm_generate_opts(model, parameters, 1), collapse = ' ')

      if (requires_files(model)) file <- tempfile('scrm')
      else file <- ''

      sum_stats <- scrm(args, file)

      seg_sites <- lapply(sum_stats$seg_sites, function(x) {
        attr(x, 'positions') <- as.numeric(colnames(x))
        x
      })

      sum_stats <- calc_sumstats(seg_sites, file, model, parameters)
      unlink(file)

      sum_stats
    },
    get_cmd = function(model) {
      cmd <- scrm_generate_opts(model, NULL, 1, FALSE)
      paste("scrm", paste(cmd, collapse = ' '))
    }
  )
)


register_simulator(SimulatorScrm) #nolint
