#' @importFrom scrm scrm
#' @importFrom phyclust ms
#' @include simulator_class.R
Simulator_scrm <- R6Class('Simulator_scrm', inherit = Simulator,
  private = list(
    name = 'scrm',
    features = c("sample", "mutation", "migration", "migration_sym",
                  "pop_merge", "recombination", "size_change", "growth"),
    sumstats = c("jsfs", "seg.sites", "file", "trees"),
    priority = 90
  ),
  public = list(
    simulate = function(dm, parameters) {
      args <- paste(sum(get_sample_size(dm, for_sim = TRUE)),
                    get_locus_number(dm),
                    paste(ms_generate_opts(dm, parameters), collapse = ' '))

      if ('file' %in% get_summary_statistics(dm)) file <- tempfile('scrm')
      else file <- ''

      sum_stats <- scrm(args, file)

      seg_sites <- lapply(sum_stats$seg_sites, function(x) {
        attr(x, 'positions') <- as.numeric(colnames(x))
        x
      })

      sum_stats <- calc_sumstats(seg_sites, file, dm, parameters)
      unlink(file)

      sum_stats
    },
    get_cmd = get_simulator("ms")$get_cmd
  )
)


register_simulator(Simulator_scrm)
