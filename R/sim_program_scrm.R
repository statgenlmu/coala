scrm_features  <- c("sample", "mutation", "migration", "migration_sym",
                    "pop_merge", "recombination", "size_change", "growth")
scrm_sum_stats <- c("jsfs", "seg.sites", "file", "trees")

#' @importFrom scrm scrm
scrm_simulate <- function(dm, parameters) {
  args <- paste(sum(get_sample_size(dm)),
                get_locus_number(dm),
                paste(generateMsOptions(dm, parameters), collapse = ' '))

  if ('file' %in% get_summary_statistics(dm)) {
    file <- tempfile('csr_scrm')
    sum_stats <- scrm(args, file)
    return(generateSumStats(file, 'ms', parameters, dm))
  }

  sum_stats <- scrm(args)
  generateSumStats(sum_stats, "scrm", parameters, dm)
}


#' @include sim_program.R
createSimProgram("scrm", scrm_features, scrm_sum_stats,
                  scrm_simulate, ms_get_command, 90)
rm(scrm_features, scrm_sum_stats, scrm_simulate)
