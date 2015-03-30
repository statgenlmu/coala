#' @importFrom scrm scrm
#' @include simulator_class.R
#' @include simulator_ms.R
SimulatorScrm <- R6Class('SimulatorScrm', inherit = Simulator, #nolint
  private = list(
    name = 'scrm',
    features = c("sample", "mutation", "migration", "migration_sym",
                  "pop_merge", "recombination", "size_change", "growth"),
    sumstats = c("jsfs", "seg.sites", "file", "trees"),
    priority = 90
  ),
  public = list(
    simulate = function(model, parameters) {
      if (length(get_locus_group_number(model)) > 1)
        stop("scrm can only simulate one group of loci at the moment.")

      args <- paste(sum(get_sample_size(model, for_sim = TRUE)),
                    get_locus_number(model),
                    paste(ms_generate_opts(model, parameters, 1),
                          collapse = ' '))

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
    get_cmd = get_simulator("ms")$get_cmd
  )
)


register_simulator(SimulatorScrm) #nolint
