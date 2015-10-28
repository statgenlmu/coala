#' @importFrom assertthat is.number
stat_omega_class <- R6Class("stat_omega", inherit = sumstat_class,
  private = list(
    req_files = TRUE,
    binary = NULL,
    min_win = NULL,
    max_win = NULL,
    grid = NULL
  ),
  public = list(
    initialize = function(name, min_win, max_win, grid, binary) {
      assert_that(is.number(min_win))
      private$min_win <- min_win
      assert_that(is.number(max_win))
      assert_that(min_win < max_win)
      private$max_win <- max_win
      assert_that(is.number(grid))
      private$grid <- grid

      if (identical(binary, "automatic")) {
        binary <- search_executable("OmegaPlus", envir_var = "OMEGAPLUS")
        if (is.null(binary)) stop("No binary for OmegaPlus found.")
      } else {
        assert_that(length(binary) == 1)
        assert_that(is.character(binary))
        assert_that(file.exists(binary))
      }
      private$binary <- binary
      super$initialize(name)
    },
    calculate = function(seg_sites, trees, files, model) {
      if (has_trios(model)) {
        stop("OmegaPlus can not be calculated from locus trios")
      }

      cur_wd <- getwd()

      op_list <- lapply(seq(along = files), function(i) {
        tmp_dir <- tempfile("omegaprime")
        dir.create(tmp_dir)
        setwd(tmp_dir)
        system2(private$binary,
                args = c("-name op",
                         "-minwin", self$get_min_win(),
                         "-maxwin", self$get_max_win(),
                         "-grid", self$get_grid(),
                         "-length", get_locus_length(model, group = i),
                         "-input", files[i]),
                stdout = TRUE)
        op <- self$parse_report(tmp_dir, self$get_grid(), 1)
        unlink(tmp_dir, recursive = TRUE)
        op
      })

      setwd(cur_wd)

      do.call(rbind, op_list)
    },
    parse_report = function(dir, n_grid, start_locus) {
      readLines(file.path(dir, "OmegaPlus_Report.op"))
      values <- read.delim(file.path(dir, "OmegaPlus_Report.op"), header = FALSE,
                           comment.char = "/")
      colnames(values) <- c("pos", "omega")
      assert_that(nrow(values) %% n_grid == 0)
      data.frame(locus = rep(start_locus:(nrow(values) / n_grid), each = n_grid),
                 values)
    },
    get_grid = function() private$grid,
    get_min_win = function() private$min_win,
    get_max_win = function() private$max_win
  )
)

#' Returns the Segregation Sites Statistics from simulations
#'
#' @inheritParams sumstat_four_gamete
#' @export
sumstat_omega <- function(name = "omega", min_win = 100, max_win = 1000,
                              grid = 10000, binary = "automatic") {
  stat_omega_class$new(name, min_win, max_win, grid, binary)
}

has_omega <- function() {
  !is.null(search_executable("OmegaPlus", envir_var = "OMEGAPLUS"))
}
