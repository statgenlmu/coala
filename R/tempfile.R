tempfile <- function(name='unnamed') {
  base::tempfile(paste0('coalsimr-', Sys.getpid(), '-', name, '-'))
}

