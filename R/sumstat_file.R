SumstatFile <- R6Class('SumstatFile', inherit = Sumstat, #nolint
  private = list(
    folder = NULL,
    req_files = TRUE
  ),
  public = list(
    initialize = function(folder) {
      dir.create(folder, showWarnings = FALSE)
      private$folder <- folder
      super$initialize('file')
    },
    calculate = function(seg_sites, files, model) {
      if (is.list(files)) files <- unlist(files)
      if (!all(file.copy(files, private$folder, overwrite = FALSE)))
        stop('Failed to copy simulated files. Look for warnings.')

      file.path(private$folder, basename(files))
    }
  )
)


#' Returns files with the raw results of simulations
#'
#' @param group The locus group for which this summary statistic is reported.
#'   The default of `0` corresponds to all groups.
#' @export
sumstat_file <- function(folder) {
  SumstatFile$new(folder) #nolint
}
