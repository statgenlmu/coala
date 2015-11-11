#' @importFrom R6 R6Class
locus_class <- R6Class("locus",
  private = list(
    number = 0,
    length = NA
  ),
  public = list(
    initialize = function(locus_length, locus_number = 1) {
      assert_that(is.numeric(locus_length))
      assert_that(length(locus_length) == 1 || length(locus_length) == 5)
      assert_that(all(locus_length >= 0))
      private$length <- locus_length

      assert_that(is.numeric(locus_number))
      assert_that(locus_number > 0)
      assert_that(length(locus_number) == 1)
      private$number <- locus_number
    },
    get_length = function(trios = FALSE) {
      if (!trios) return(sum(private$length))

      if (length(private$length) == 1) {
        locus_length <- c(0, 0, private$length, 0, 0)
      } else if (length(private$length) == 5) {
        locus_length <- private$length
      } else stop("Failed to get locus length")

      names(locus_length) <- c("length_l", "length_il", "length_m",
                               "length_ir", "length_r")
      locus_length
    },
    get_number = function() private$number,
    print = function() {
      cat(self$get_number(), ifelse(self$get_number() == 1, "locus", "loci"),
          "of length", self$get_length(), "\n")
    }
  )
)


is.locus <- function(locus) any("locus" == class(locus))


#' Add one locus or multiple loci to a Model
#'
#' @param length The length of the locus in base pairs.
#' @export
#' @examples
#' # A model with one locus of length 1005 bp
#' coal_model(5:7, 0) + locus_single(1005)
#'
#' # A model with ten loci of average length 950bp
#' coal_model(15, 0) + locus_averaged(10, 950)
#' # or just
#' coal_model(15, 10, 950)
#'
#' # A model with two loci. The first group consists of 10 loci with
#' # a length 560bp each, the second one of two loci with length 750bp and 560pb,
#' # respectively.
#' coal_model(20, 10, 560) +
#'   locus_single(750) +
#'   locus_single(430)
locus_single <- function(length) {
  locus_class$new(length, 1)
}


#' @describeIn locus_single Multiple Loci of the same length.
#' @param number The number of loci to add.
#' @export
locus_averaged <- function(number, length) {
  locus_class$new(round(length), number)
}


#' Adds a trio of loci to a group
#'
#' @inheritParams locus_single
#' @param locus_length An integer vector of length 3, giving the length of each
#'   of the three loci (left, middle and right).
#' @param distance A vector of two, giving the distance between left and middle,
#'   and middle an right locus, in base pairs.
#' @export
#' @examples
#' coal_model(c(25,25)) +
#'   locus_trio(locus_length=c(1250, 1017, 980),
#'              distance=c(257, 814))
locus_trio <- function(locus_length = c(left = 1000,
                                        middle = 1000,
                                        right = 1000),
                       distance = c(left_middle = 500,
                                    middle_right = 500),
                       number = 1) {

  stopifnot(length(locus_length) == 3)
  if (!is.null(names(locus_length))) {
    locus_length <- locus_length[c("left", "middle", "right")]
  }
  stopifnot(length(distance) == 2)
  if (!is.null(names(distance))) {
    distance <- distance[c("left_middle", "middle_right")]
  }

  locus_class$new(locus_length = c(locus_length, distance)[c(1, 4, 2, 5, 3)],
                  locus_number = number)
}


# Converts a position on the middle locus to the relative position
# on the simulated stretch
conv_middle_to_trio_pos <- function(pos, model, locus,
                                    relative_out = TRUE, relative_in = TRUE) {

  llm <- get_locus_length_matrix(model)
  group <- get_locus_group(model, locus)

  if (relative_in) pos <- pos * llm[group, 3]
  pos <- pos + llm[group, 1] + llm[group, 2]
  if (relative_out) pos <- pos / sum(llm[group, 1:5])

  names(pos) <- NULL
  pos
}
