#' @importFrom R6 R6Class
SumstatNsl <- R6Class('sumstat_nsl', inherit = SumstatIhh, #nolint
  public = list(
    segsites_to_snp_map = function(seg_sites, pos) {
      map <- data.frame(name = seq(along = pos),
                        chr = 1,
                        pos = seq(along = pos),
                        anc = 0,
                        der = 1)
      file <- tempfile('snp_map')
      write.table(map, file, row.names = FALSE, col.names = FALSE)
      file
    },
    calculate = function(seg_sites, files, model) {
      stat <- super$calculate(seg_sites, files, model)
      lapply(stat, function(x) {
        colnames(x) <- NULL
        x[ , 3]
      })
    }
  )
)

#' @export
sumstat_nSL <- function(name = 'ihh', position=NA, population=1) { #nolint
  SumstatNsl$new(name, position, population) #nolint
}
