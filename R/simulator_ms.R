ms_features  <- c("sample", "mutation", "migration", "migration_sym",
                  "pop_merge", "recombination", "size_change", "growth",
                  "inter_locus_variation", "trees",
                  "unphased", "ploidy", "samples_per_ind")
ms_sum_stats <- c("jsfs", "trees", "seg.sites", "file")


# This function generates an string that contains an R command for generating
# an ms call to the current model.
ms_generate_opts_cmd <- function(dm) {
  sample_size <- get_sample_size(dm, for_sim = TRUE)
  cmd <- c('c(')
  cmd <- c(cmd,'"-I"', ",", length(sample_size), ',',
           paste(sample_size, collapse=","), ',')

  for (i in 1:dim(get_feature_table(dm))[1] ) {
    type <- as.character(get_feature_table(dm)[i,"type"])
    feat <- unlist(get_feature_table(dm)[i, ])

    if ( type == "mutation" ) {
      cmd <- c(cmd,'"-t"', ',', feat["parameter"], ',')
    }

    else if (type == "pop_merge") {
      cmd <- c(cmd, '"-ej"', ',', feat["time.point"], ',',
               feat["pop.source"], ',', feat["pop.sink"], ',')
    }

    else if (type == "migration")
      cmd <- c(cmd, '"-em"', ',', feat['time.point'], ',',
               feat['pop.sink'], ',', feat['pop.source']  , ',',
               feat['parameter'], ',')

    else if (type == "migration_sym")
      cmd <- c(cmd, '"-eM"', ',',
               feat['time.point'], ',',
               feat['parameter'], ',')

    else if (type == "recombination")
      cmd <- c(cmd, '"-r"', ',', feat['parameter'], ',',
               get_locus_length(dm), ',')

    else if (type == "size_change"){
      cmd <- c(cmd, '"-en"', ',', feat['time.point'], ',',
               feat["pop.source"], ',', feat['parameter'], ',')
    }

    else if (type == "growth"){
      cmd <- c(cmd, '"-eg"', ',' , feat["time.point"], ',',
               feat["pop.source"], ',', feat["parameter"], ',')
      }

    else if (type == 'trees') {
      cmd <- c(cmd, '"-T",')
    }

    else if (type %in% c("sample", "loci.number", "loci.length",
                         "selection", "selection_AA", "selection_Aa",
                         "inter_locus_variation", "unphased",
                         "ploidy", "samples_per_ind")) NULL
    else stop("Unknown feature:", type)
  }


  cmd <- c(cmd, '" ")')
}

ms_generate_opts <- function(dm, parameters, eval_pars = TRUE) {
  ms.tmp <- create_par_env(dm, parameters)

  cmd <- read_cache(dm, 'ms_cmd')
  if (is.null(cmd)) {
    cmd <- ms_generate_opts_cmd(dm)
    cache(dm, 'ms_cmd', cmd)
  }

  if (!eval_pars) cmd <- escape_par_expr(cmd)
  eval(parse(text=cmd), envir=ms.tmp)
}


#' @importFrom phyclust ms
#' @include simulator_class.R
Simulator_ms <- R6Class('Simulator_ms', inherit = Simulator,
  private = list(
    name = 'ms',
    features = ms_features,
    sumstats = ms_sum_stats,
    priority = 100
  ),
  public = list(
    simulate = function(dm, parameters=numeric()) {
      stopifnot(length(parameters) == 0 | all(is.numeric(parameters)))

      # Run all simulation in with one ms call if they loci are identical,
      # or call ms for each locus if there is variation between the loci.
      if (has_inter_locus_var(dm)) {
        sim_reps <- 1:get_locus_number(dm)
        sim_loci <- 1
      } else {
        sim_reps <- 1
        sim_loci <- get_locus_number(dm)
      }

      # Do the actuall simulation
      files <- lapply(sim_reps, function(locus) {
        ms.options <- ms_generate_opts(dm, parameters, locus)
        file <- tempfile('csr_ms')

        ms(sum(get_sample_size(dm, for_sim = TRUE)), sim_loci,
           unlist(strsplit(ms.options, " ")), file)

        if(file.info(file)$size == 0) stop("ms simulation output is empty")
        file
      })

      # Parse the output and calculate summary statistics
      seg_sites <- parse_ms_output(files,
                                   get_sample_size(dm, for_sim = TRUE),
                                   get_locus_number(dm))

      sum_stats <- calc_sumstats(seg_sites, files, dm, parameters)

      # Clean Up
      unlink(files)
      sum_stats
    },
    get_cmd = function(model) {
      cmd <- ms_generate_opts(model, get_parameter_table(model)$name, FALSE)
      txt <- paste(cmd, collapse = ' ')
      paste("ms", sum(get_sample_size(model)), get_locus_number(model), txt)
    }
  )
)


register_simulator(Simulator_ms)
