sg_mutation_models <- c("HKY", "GTR")


generate_tree_model <- function(model) {
  tree_model <- read_cache(model, "tree_model")

  if (is.null(tree_model)) {
    tree_model <- model

    # Features
    tree_model_features <- !vapply(model$features, function(x) {
      any(c("seg_sites_feat",
            "mutation",
            "outgroup") %in% class(x))
    }, logical(1)) #nolint
    if (all(tree_model_features)) stop("seq-gen not required")
    tree_model$features <- model$features[tree_model_features]

    # Summary Stastics
    tree_model$sum_stats <- create_sumstat_container()
    tree_model <- tree_model + sumstat_trees()

    cache(model, "tree_model", tree_model)
  }

  tree_model
}


conv_to_seqgen_arg <- function(feature, model) UseMethod("conv_to_seqgen_arg")

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_seqgen_arg.default <- function(feature, model) {
  stop("Unknown feature when generating seqgen command")
}


sg_generate_opts <- function(model, parameters, locus,
                             seed = NULL, for_cmd = FALSE) {
  locus_lengths <- get_locus_length(model, group = locus, total = FALSE)

  if (length(locus_lengths) == 5) {
    locus_lengths <- locus_lengths[c(1, 3, 5)]
  }

  cmd <- get_simulator("seqgen")$create_cmd_template(model)
  #print(cmd)

  # Fill the parameters in the template
  vapply(seq_along(locus_lengths), function(i) {
    par_envir <- create_par_env(model, parameters, locus = locus,
                                locus_length = locus_lengths[i],
                                for_cmd = for_cmd)
    paste(eval(parse(text = cmd[[i]]), envir = par_envir), collapse = " ")
  }, character(1L))
}


#' @importFrom R6 R6Class
#' @include simulator_class.R
seqgen_class <- R6Class("seqgen", inherit = simulator_class,
  private = list(
    name = "seqgen",
    binary = NULL,
    priority = 100
  ),
  public = list(
    initialize = function(binary = NULL, priority = 100) {
      # Try to automatically find a jar file and java if not given
      if (is.null(binary)) {
        binary <- search_executable(c("seqgen", "seq-gen",
                                      "seqgen.exe", "seq-gen.exe"), "SEQGEN")
      }
      if (is.null(binary)) stop("No binary file for seqgen found.")
      if (!file.exists(binary)) stop("seqgen binary (", binary,
                                     ") does not exist.")
      message("Using '", binary, "' as seqgen binary")
      assert_that(is.character(binary) && length(binary) == 1)
      private$binary <- binary

      super$initialize(priority)
    },
    create_task = function(model, pars, locus_number,
                           locus_id = 1,
                           eval_pars = TRUE) {
      tree_model <- generate_tree_model(model)
      tree_simulator <- select_simprog(tree_model)
      tree_task <- tree_simulator$create_task(tree_model, pars, locus_number,
                                              locus_id, eval_pars)

      cmd <- sg_generate_opts(model, pars, locus_id)
      trio_dists <- get_locus_length_matrix(model)[1, 1:5]

      create_sim_task(self, locus_number,
                      cmd = cmd,
                      trio_dists = trio_dists,
                      tree_model = tree_model,
                      tree_task = tree_task)
    },
    create_cmd_template = function(model) {
      cmd <- read_cache(model, "seqgen_cmd")

      if (is.null(cmd)) {
        if (has_trios(model)) is_outer <- c(TRUE, FALSE, TRUE)
        else is_outer <- FALSE

        cmd <- lapply(is_outer, function(outer) {
          cmd <- paste(vapply(model$features, conv_to_seqgen_arg,
                              FUN.VALUE = character(1), model = model),
                       collapse = "")
          cmd <- paste0("c('", cmd, "')")
        })

        cache(model, "seqgen_cmd", cmd)
      }
      cmd
    },
    call = function(args) {
      suppressWarnings(results <- system2(private$binary,
                                          paste(args, "-z", sample_seed(1)),
                                          stdout = TRUE))
      results
    },
    simulate = function(model, sim_task) {
      assert_that(is.model(model))
      assert_that(is_simulation_task(sim_task))

      trio_dists <- sim_task$get_arg("trio_dists")
      has_trios <- sum(trio_dists > 0) > 1

      # Simulate the ancestral trees
      tree_task <- sim_task$get_arg("tree_task")
      tree_task$locus_number <- sim_task$locus_number
      tree_model <- sim_task$get_arg("tree_model")
      trees <- tree_task$get_simulator()$simulate(tree_model, tree_task)$trees
      assert_that(is.list(trees))
      assert_that(length(trees) == tree_task$locus_number)

      # Prepare the simulation commands
      cmds <- sim_task$get_arg("cmd")

      # Prepare the trees for seqgen
      tree_files <- tempfile(c(left = "left",
                               middle = "middle",
                               right = "right"))
      generate_trio_trees(trees, trio_dists, tree_files)

      if (!has_trios) {
        unlink(tree_files[-2])
        tree_files <- tree_files[2]
        locus_length <- trio_dists[3]
      } else {
        locus_length <- trio_dists[c(1, 3, 5)]
      }

      assert_that(length(tree_files) == length(cmds))
      assert_that(length(tree_files) == length(locus_length))

      # Call seq-gen for each trio locus
      sim_results <- lapply(seq_along(cmds), function(trio_locus) {
        cmd <- paste(cmds[trio_locus], tree_files[trio_locus])

        sim_output <- self$call(cmd)
        result <- list(cmd = cmd)

        if (requires_segsites(model)) {
          result$seg_sites <-
            parse_seqgen_output(sim_output,
                                individuals = sum(get_sample_size(model, TRUE)),
                                locus_length = locus_length[trio_locus],
                                locus_number = sim_task$locus_number,
                                outgroup_size = get_outgroup_size(model, TRUE),
                                calc_segsites = TRUE)
        }

        if (requires_files(model)) {
          result$files <- tempfile("seqgen_result")
          write(sim_output, file = result$files, sep = "\n")
        }

        result
      })
      unlink(tree_files)

      if (length(sim_results) == 1) {
        output <- sim_results[[1]]
      } else {
        assert_that(length(sim_results) == 3)
        output <- list(
          seg_sites = create_locus_trio(sim_results[[1]]$seg_sites,
                                        sim_results[[2]]$seg_sites,
                                        sim_results[[3]]$seg_sites),
          files = c(sim_results[[1]]$files,
                    sim_results[[2]]$files,
                    sim_results[[3]]$files)
        )
      }

      output$simulator <- self
      output
    },
    get_cmd = function(model) {
      c(trees = get_cmd(generate_tree_model(model)),
        sequence = paste("seqgen", sg_generate_opts(model, NULL, 1, 0, TRUE),
                         collapse = " "))
    },
    get_info = function() c(name = "seqgen", binary = private$binary)
  )
)

has_seqgen <- function() !is.null(simulators[["seqgen"]])


#' Simulator: seq-gen
#'
#' This allows you to use seq-gen to simulate finite sites mutation models.
#' When using seq-gen, coala will simulate ancestral tress using the other
#' simulators, and call seq-gen to simulate finite sites mutations using the
#' trees. Seq-gen has a low priority, but will always be used when finite
#' sites mutation models are used.
#'
#' @section Installation:
#' You need to download the program from
#' \url{http://tree.bio.ed.ac.uk/software/seqgen/}
#' and compile the binary prior to invoking this function.
#' On Debian-based systems, you can alternatively install the package
#' 'seg-gen'.
#'
#' @references
#' Andrew Rambaut and Nicholas C. Grassly.
#' Seq-Gen: an application for the Monte Carlo simulation of DNA sequence
#' evolution along phylogenetic trees.
#' Comput Appl Biosci (1997) 13 (3): 235-238
#' doi:10.1093/bioinformatics/13.3.235
#'
#' @param binary The path of the seqgen binary that will be used
#'  for simulations. If none is provided, coala will look for a
#'  binary called 'seqgen' or 'seq-gen' using the PATH variable.
#' @inheritParams simulator_ms
#' @name simulator_seqgen
#' @family simulators
#' @export
#' @examples
#' \dontrun{activate_seqgen("./bin/seqgen")}
activate_seqgen <- function(binary = NULL, priority = 100) {
  register_simulator(seqgen_class$new(binary, priority))
  reset_cache()
  invisible(NULL)
}
