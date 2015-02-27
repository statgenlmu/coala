#' @importFrom R6 R6Class
Locus <- R6Class('Locus', inherit = base_class,
  private = list(
    group = 0,
    length = NA,
    number = NA,
    name = NA
  ),
  public = list(
    initialize = function(length, number = 1, name = '', group = 0) {
      stopifnot(is.numeric(length))
      stopifnot(length > 0)
      private$length <- length
      stopifnot(is.numeric(number))
      stopifnot(number > 0)
      private$number <- number
      stopifnot(is.character(name))
      private$name <- name
      stopifnot(is.numeric(group))
      stopifnot(group >= 0)
      private$group <- group
    },
    get_name = function() private$name,
    get_length = function() private$length,
    get_number = function() private$number,
    get_group = function() private$group
  )
)


is.locus <- function(locus) 'Locus' %in% class(locus)


#' Add one locus or multiple loci to a Model
#'
#' @param length The length of the locus in base pairs.
#' @param group The group of loci to which the locus belongs.
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
#' # A model with two groups of loci. The first group consists of 10 loci with
#' # a length 560bp each, the second one of two loci with length 750bp and 560pb,
#' # repectively.
#' coal_model(20, 10, 560) +
#'   locus_single(750, group = 2) +
#'   locus_single(430, group = 2)
locus_single <- function(length, group = 0, name = '') {
  Locus$new(length, 1, name = name, group = group)
}


#' @describeIn locus_single
#' @param number The number of loci to add.
#' @export
locus_averaged <- function(number, length, group = 0) {
  Locus$new(round(length), number, group = group)
}


#' Adds a trio of loci to a group
#'
#' @inheritParams locus_single
#' @param locus_names A vector of 3 strings, giving the names for the loci.
#'   The names are used for identifying the loci later (left, middle and right).
#' @param locus_length An integer vector of length 3, giving the length of each
#'   of the three loci (left, middle and right).
#' @param distance A vector of two, giving the distance between left and middle,
#'   and middle an right locus, in basepairs.
#' @export
#' @examples
#' coal_model(c(25,25)) +
#'   locus_trio(locus_names = c('Solyc00g00500.2',
#'                              'Solyc00g00520.1',
#'                              'Solyc00g00540.1'),
#'              locus_length=c(1250, 1017, 980),
#'              distance=c(257, 814))
locus_trio <- function(locus_names = c(left='', middle='', right=''),
                       locus_length = c(left=1000, middle=1000, right=1000),
                       distance = c(left_middle=500, middle_right=500),
                       number = 1,
                       group = 0) {

  stopifnot(length(locus_names) == 3)
  if (!is.null(names(locus_names))) {
    locus_names <- locus_names[c('left', 'middle', 'right')]
  }
  stopifnot(length(locus_length) == 3)
  if (!is.null(names(locus_length))) {
    locus_length <- locus_length[c('left', 'middle', 'right')]
  }
  stopifnot(length(distance) == 2)
  if (!is.null(names(distance))) {
    distance <- distance[c('left_middle', 'middle_right')]
  }

  Locus$new(length = c(locus_length, distance)[c(1,4,2,5,3)],
            number = number,
            name = locus_names,
            group = group)
}
