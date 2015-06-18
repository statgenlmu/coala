sg_mutation_models <- c('HKY', 'GTR')


generate_tree_model <- function(model) {
  tree_model <- read_cache(model, "tree_model")

  if (is.null(tree_model)) {
    tree_model <- model

    # Features
    tree_model_features <- !vapply(model$features, function(x) {
      any(c("Feature_seg_sites",
            "Feature_mutation",
            "Feature_outgroup") %in% class(x))
    }, logical(1)) #nolint
    if (all(tree_model_features)) stop("seq-gen not required")
    tree_model$features <- model$features[tree_model_features]

    # Summary Stastics
    tree_model$sum_stats <- create_sumstat_container()
    tree_model <- tree_model + sumstat_sg_trees()

    cache(model, "tree_model", tree_model)
  }

  tree_model
}


# Function to perform simulation using seqgen
#
# @param opts The options to pass to ms. Must either be a character or character
# vector.
sg_call <- function(opts, ms_files) {
  stopifnot(length(opts) == length(ms_files))

  sapply(seq(along = opts), function(i) {
    if (file.info(ms_files[[i]])$size <= 1 ) {
      stop("No trees in file ", ms_files[[i]])
    }

    seqgen_file <- tempfile('seqgen')
    cmd <- paste(opts[i], ms_files[[i]])

    ret <- system2(get_executable("seqgen"), cmd, stdout = seqgen_file)

    if (ret != 0) stop("seq-gen simulation failed")
    if (!file.exists(seqgen_file)) stop("seq-gen simulation failed!")
    if (file.info(seqgen_file)$size <= 1) stop("seq-gen output is empty!")

    seqgen_file
  })
}


conv_to_seqgen_arg <- function(feature, model) UseMethod("conv_to_seqgen_arg")

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_seqgen_arg.default <- function(feature, model) {
  stop("Unknown feature when generating seqgen command")
}


sg_generate_opts <- function(model, parameters, locus,
                             seeds, for_cmd = FALSE) {
  locus_lengths <- get_locus_length(model, group = locus, total = FALSE)

  if (length(locus_lengths) == 5) {
    locus_lengths <- locus_lengths[c(1, 3, 5)]
  }

  cmd <- sg_generate_opt_cmd(model)
  #print(cmd)

  # Fill the parameters in the template
  sapply(seq(along = locus_lengths), function(i) {
    par_envir <- create_par_env(model, parameters, locus = locus,
                                locus_length = locus_lengths[i],
                                seed = seeds[i], for_cmd = for_cmd)
    paste(eval(parse(text = cmd[[i]]), envir = par_envir), collapse = " ")
  })
}


sg_generate_opt_cmd <- function(model) {

  cmd <- read_cache(model, 'seqgen_cmd')

  if (is.null(cmd)) {
    if (has_trios(model)) is_outer <- c(TRUE, FALSE, TRUE)
    else is_outer <- FALSE

    cmd <- lapply(is_outer, function(outer) {
      cmd <- paste(vapply(model$features, conv_to_seqgen_arg,
                          FUN.VALUE = character(1), model = model),
                   collapse = "")
      cmd <- paste0("c('", cmd, "')")
    })

    cache(model, 'seqgen_cmd', cmd)
  }
  cmd
}


sg_simulate <- function(model, parameters) {
  # Simulate the ancestral trees
  tree_model <- generate_tree_model(model)
  trees <- simulate(tree_model, pars = parameters)$trees
  assert_that(!is.null(trees))

  # Call seq-gen for each locus (trio)
  seqgen_files <- lapply(1:length(trees), function(locus) {
    seqgen_options <- sg_generate_opts(model, parameters, locus,
                                       sample_seed(length(trees[[locus]])))
    sg_call(seqgen_options, trees[[locus]])
  })
  unlink(unlist(trees))

  # Generate the summary statistics
  if (requires_segsites(model)) {
    seg_sites <- parse_sg_output(seqgen_files,
                                 sum(get_sample_size(model, for_sim = TRUE)),
                                 get_locus_length_matrix(model),
                                 get_locus_number(model),
                                 outgroup_size = get_outgroup_size(model, TRUE))
  } else {
    seg_sites <- NULL
  }

  sum_stats <- calc_sumstats(seg_sites, seqgen_files, model, parameters, NULL,
                             get_simulator("seqgen"))

  # Clean Up
  unlink(unlist(seqgen_files))
  sum_stats
}


#' @importFrom R6 R6Class
#' @include simulator_class.R
Simulator_seqgen <- R6Class('Simulator_seqgen', inherit = Simulator,
  private = list(
    name = 'seqgen',
    priority = 10
  ),
  public = list(
    simulate = sg_simulate,
    get_cmd = function(model) {
      c(trees = get_cmd(generate_tree_model(model)),
        sequence = paste("seqgen", sg_generate_opts(model, NULL, 1, 0, TRUE),
                         collapse = ' '))
    }
  )
)


register_simulator(Simulator_seqgen)
