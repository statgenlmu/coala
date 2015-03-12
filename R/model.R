#' @export
coal_model <- function(sample_size=0, loci_number=0, loci_length=1000) {
  model <- list()
  class(model) <- c("Coalmodel")

  model$features <- list()
  model$loci <- list()
  model$parameter <- list()
  model$sum_stats <- create_sumstat_container()

  model$scaling_factor <- 1
  model$id <- get_id()

  # Add sample sizes
  for (pop in seq(along = sample_size)) {
    if (sample_size[pop] > 0) model <- model +
      feat_sample(sample_size[pop], pop)
  }

  # Add locus
  if (loci_number > 0) {
    model <- model + locus_averaged(loci_number, loci_length)
  }

  model
}


is.model <- function(model) {
  "Coalmodel" %in% class(model)
}


# Checks if a vector of parameters is within the ranges of the model
check_par_range <- function(dm, param) {
  if (length(param) != nrow(get_parameter_table(dm))) {
    stop("Wrong number of parameters")
  }

  ranges <- get_parameter_table(dm)[,2:3]
  in.range <- all(ranges[, 1] - 1e-11 <= param & param <= ranges[, 2] + 1e-11)
  if (!in.range) stop("Parameter combination out of range")
}

# Selects a program for simulation that is capable of all current features
determine_simprog <- function(dm) {
  name <- read_cache(dm, 'simprog')

  if (is.null(name)) {
    priority <- -Inf

    for (simprog_name in ls(simulators)) {
      simprog <- get_simulator(simprog_name)
      if (all(get_feature_table(dm)$type %in% simprog$get_features()) &
            all(dm$sum_stats$name %in% simprog$get_sumstats())) {

        if (simprog$get_priority() > priority) {
          name <- simprog
          priority <- simprog$get_priority()
        }
      }
    }

    if (is.null(name)) stop("No suitable simulation software found!")
    cache(dm, 'simprog', name)
  }

  name
}


get_mutation_par <- function(dm, outer=FALSE) {
  if (outer) {
    feat <- search_feature(dm, "mutation_outer")
    if (nrow(feat) == 0) {
      feat <- search_feature(dm, "mutation")
    }
  }  else {
    feat <- search_feature(dm, "mutation")
  }
  if (nrow(feat) != 1) stop("Failed to determine mutation rate")
  feat[1, 'parameter']
}


search_feature <- function(dm, type=NULL, pop.source=NULL,
                          pop.sink=NULL, time.point=NULL,
                          feat_table=TRUE) {

  feat_tbl <- get_feature_table(dm)
  mask <- rep(TRUE, nrow(feat_tbl))

  if (!is.null(type)) mask <- mask & feat_tbl$type %in% type

  if (!is.null(pop.source)) {
    if (is.na(pop.source)) {
      mask <- mask & is.na(feat_tbl$pop.source)
    } else {
      mask <- mask & feat_tbl$pop.source %in% pop.source
    }
  }

  if (!is.null(pop.sink)) {
    if (is.na(pop.sink)) {
      mask <- mask & is.na(feat_tbl$pop.sink)
    } else {
      mask <- mask & feat_tbl$pop.sink %in% pop.sink
    }
  }

  if (!is.null(time.point)) {
    if (is.na(time.point)) {
      mask <- mask & is.na(feat_tbl$time.point)
    } else {
      mask <- mask & feat_tbl$time.point %in% time.point
    }
  }

  if(!feat_table) return(get_features(dm)[mask])
  feat_tbl[mask, ]
}


add_inter_locus_var <- function(dm) {
  if (has_inter_locus_var(dm)) return(dm)
  dm + Feature$new('inter_locus_variation', par_const(NA))
}


has_inter_locus_var <- function(dm) {
  nrow(search_feature(dm, 'inter_locus_variation')) > 0
}


has_trios <- function(dm) {
  sum(get_locus_length_matrix(dm)[,-3]) > 0
}


get_snp_positions <- function(seg_sites, model, relative=TRUE) {
  assert_that(length(seg_sites) == get_locus_number(model))
  llm <- get_locus_length_matrix(model)
  lapply(1:length(seg_sites), function(locus) {
    pos <- attr(seg_sites[[locus]], 'position')
    trio_locus <- attr(seg_sites[[locus]], 'locus')
    if (is.null(trio_locus)) trio_locus <- 0
    pos[trio_locus == -1] <- pos[trio_locus == -1] * llm[locus, 1]
    pos[trio_locus == 0] <- pos[trio_locus == 0] * llm[locus, 3] +
      sum(llm[locus, 1:2])
    pos[trio_locus == 1] <- pos[trio_locus == 1] * llm[locus, 5] +
      sum(llm[locus, 1:4])
    if (relative) pos <- pos / sum(llm[locus, 1:5])
    pos
  })
}
