SumstatTrees <- R6Class('SumstatTrees', inherit = Sumstat, #nolint
  private = list(req_files = TRUE),
  public = list(
    calculate = function(seg_sites, files, model) {
      parse_trees(files, get_locus_number(model))
    }
  )
)

#' Returns ancestral tress in NEWICK fromat from simulations
#'
#' @export
#' @inheritParams sumstat_file
sumstat_trees <- function(name = "trees") {
  coal_model() +
    Feature$new('trees', par_const(NA)) +
    SumstatTrees$new(name) #nolint
}


SumstatSgTrees <- R6Class('SumstatSgTrees', inherit = Sumstat, #nolint
  private = list(req_files = TRUE),
  public = list(
    calculate = function(seg_sites, files, model) {
      trees <- parse_trees(files, get_locus_number(model), FALSE)
      llm <- get_locus_length_matrix(model)

      trio_trees <- generate_trio_trees(trees, llm)
      lapply(trio_trees, function(locus_trees) {
        files <- sapply(locus_trees, function(locus) {
          if (length(locus) > 0) {
            file <- tempfile("trio_tree")
            write(locus, file, sep="\n")
          } else {
            file <- NULL
          }
          file
        })
        files[!sapply(files, is.null)]
      })
    }
  )
)

#' Returns ancestral tress as files for seq-gen
#'
#' @inheritParams sumstat_file
sumstat_sg_trees <- function() {
  coal_model() +
    Feature$new('trees', par_const(NA)) +
    SumstatSgTrees$new('trees') #nolint
}
