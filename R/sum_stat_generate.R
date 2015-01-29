generateSumStats <- function(files, program, parameters, dm, seg_sites) {
  model_stats <- get_summary_statistics(dm)

  calc_seg_sites <- any(c('seg.sites', 'jsfs') %in% model_stats)
  if (missing(seg_sites) & calc_seg_sites) {
    if (program == 'ms') {
      seg_sites <- parseMsOutput(files,
                                 get_sample_size(dm),
                                 get_locus_number(dm))
    } else if (program == 'seqgen') {
      seg_sites <- parseSeqgenOutput(files,
                                     sum(get_sample_size(dm)),
                                     get_locus_length_matrix(dm),
                                     get_locus_number(dm),
                                     outgroup_size = get_outgroup_size(dm))
    } else if (program == 'scrm') {
      seg_sites <- files[['seg_sites']]
      seg_sites <- lapply(seg_sites, function(x) {
        attr(x, 'positions') <- as.numeric(colnames(x))
        x
      })
    } else {
      stop("Unknown program: ", program)
    }
  }

  # Add the parameters of the simulation
  sum_stats <- list()
  if (!missing(parameters)) sum_stats[['pars']] <- parameters

  for (sum_stat in get_summary_statistics(dm)) {
    # Add seg_sites
    if (sum_stat == 'seg.sites') {
      sum_stats[['seg.sites']] <- seg_sites
    }

    # Add JSFS
    else if (sum_stat == 'jsfs') {
      sum_stats[['jsfs']] <- calcJsfs(seg_sites, get_sample_size(dm))
    }

    else if (sum_stat == 'file' | sum_stat == 'trees') { }
    else stop('Unknown summary statistik:', sum_stat)
  }

  # Add files if needed, or delete otherwise
  if (!missing(files)) {
    if ('file' %in% model_stats) {
      sum_stats[['file']] <- files
    } else {
      unlink(unlist(files))
    }
  }

  sum_stats
}
