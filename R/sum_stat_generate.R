generateSumStats <- function(files, program, parameters, dm, seg_sites) {
  model_stats <- dm.getSummaryStatistics(dm)
  
  calc_seg_sites <- any(c('seg.sites', 'jsfs', 'pmc', 'fpc') %in% model_stats)
  if (missing(seg_sites) & calc_seg_sites) {
    if (program == 'ms') {
      seg_sites <- parseMsOutput(files, 
                                 dm.getSampleSize(dm), 
                                 dm.getLociNumber(dm))
    } else if (program == 'seqgen') {
      seg_sites <- parseSeqgenOutput(files,
                                     sum(dm.getSampleSize(dm)),
                                     dm.getLociLengthMatrix(dm),
                                     dm.getLociNumber(dm),
                                     outgroup_size = dm.getOutgroupSize(dm))
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

  for (i in 1:nrow(dm@sum.stats)) {
    stat <- dm@sum.stats[i, ]
    
    # Add seg_sites
    if (stat$name == 'seg.sites') {
      sum_stats[['seg.sites']] <- seg_sites
    }
    
    # Add JSFS
    else if (stat$name == 'jsfs') {
      sum_stats[['jsfs']] <- calcJsfs(seg_sites, dm.getSampleSize(dm))
    }
      
    else if (stat$name == 'file' | stat$name == 'trees') { }
    else stop('Unknown summary statistik:', stat$name)
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
