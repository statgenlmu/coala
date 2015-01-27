# --------------------------------------------------------------
# Uses scrm to simulate demographic models
# 
# Authors:  Lisha Mathew & Paul R. Staab
# Licence:  GPLv3 or later
# --------------------------------------------------------------

scrm_features  <- c("sample", "loci.number", "loci.length",
                    "mutation", "migration", "split",
                    "recombination", "size.change", "growth")

scrm_sum_stats <- c("jsfs", "seg.sites", "file", "trees")

scrm_simulate <- function(dm, parameters) {
  checkType(dm, "dm")
  checkType(parameters, "num")

  if (length(parameters) != dm.getNPar(dm)) stop("Wrong number of parameters!")

  args <- paste(sum(dm.getSampleSize(dm)),
                dm.getLociNumber(dm),
                paste(generateMsOptions(dm, parameters), collapse = ' '))
  
  if ('file' %in% dm.getSummaryStatistics(dm)) {
    file <- getTempFile('scrm')
    sum_stats <- scrm(args, file)
    return(generateSumStats(file, 'ms', parameters, dm))
  }

  sum_stats <- scrm(args)
  generateSumStats(sum_stats, "scrm", parameters, dm)
}

scrm_finalize <- function(dm) {
  dm@options[['ms.cmd']] <- generateMsOptionsCommand(dm)
  return(dm)
}

#' @include dm_sim_program.R
createSimProgram("scrm", scrm_features, scrm_sum_stats,
                  scrm_simulate, scrm_finalize, printMsCommand, 90)
rm(scrm_features, scrm_sum_stats, scrm_simulate, scrm_finalize)