scrm_features  <- c("sample", "mutation", "migration", "migration_sym",
                    "pop_merge", "recombination", "size_change", "growth")
scrm_sum_stats <- c("jsfs", "seg.sites", "file", "trees")

#' @importFrom scrm scrm
scrm_simulate <- function(dm, parameters) {
  args <- paste(sum(get_sample_size(dm, for_sim = TRUE)),
                get_locus_number(dm),
                paste(generateMsOptions(dm, parameters), collapse = ' '))

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
}


#' @include sim_program.R
createSimProgram("scrm", scrm_features, scrm_sum_stats,
                  scrm_simulate, ms_get_command, 90)
