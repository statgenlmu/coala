#' @importFrom R6 R6Class
SumstatSegSites <- R6Class('SumstatSegSites', inherit = Sumstat, #nolint
  private = list(req_segsites = TRUE),
  public = list(
    calculate = function(seg_sites, files, model) seg_sites
  )
)


#' Returns the Segregation Sites Statistics from simulations
#'
#' @inheritParams sumstat_file
#' @export
sumstat_seg_sites <- function(name = 'seg_sites') {
  SumstatSegSites$new(name) #nolint
}


conv_for_trios <- function(seg_sites, llm) {
  assert_that(length(seg_sites) == nrow(llm))
  for (i in seq(along = seg_sites)) {
    total_length <- sum(llm[i, 1:5])
    borders <- cumsum(llm[i, 1:4] / total_length)
    pos <- attr(seg_sites[[i]], "positions")
    left <- pos < borders[1]
    middle <- pos >= borders[2] & pos < borders[3]
    right <- pos >= borders[4]
    seg_sites[[i]] <- seg_sites[[i]][ , left | middle | right, drop=FALSE]

    pos[left] <- pos[left] * total_length / llm[i, 1]
    pos[middle] <- (pos[middle] - borders[2]) * total_length / llm[i, 3]
    pos[right] <- (pos[right] - borders[4]) * total_length / llm[i, 5]

    attr(seg_sites[[i]], "positions") <- pos[left | middle | right]
    attr(seg_sites[[i]], "locus") <- c(rep(-1, sum(left)),
                                       rep(0, sum(middle)),
                                       rep(1, sum(right)))
  }
  seg_sites
}
