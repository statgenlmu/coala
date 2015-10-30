stat_file_class <- R6Class("stat_file", inherit = sumstat_class, #nolint
  private = list(
    folder = NULL,
    req_files = TRUE
  ),
  public = list(
    initialize = function(folder) {
      dir.create(folder, showWarnings = FALSE)
      private$folder <- folder
      super$initialize("file", identity)
    },
    calculate = function(seg_sites, trees, files, model) {
      if (is.list(files)) files <- unlist(files)
      if (!all(file.copy(files, private$folder, overwrite = FALSE)))
        stop("Failed to copy simulated files. Look for warnings.")

      file.path(private$folder, basename(files))
    }
  )
)


#' Returns files with the raw results of simulations
#'
#' @param folder The path of a folder where the files will be written.
#' @export
sumstat_file <- function(folder) {
  stat_file_class$new(folder)
}
