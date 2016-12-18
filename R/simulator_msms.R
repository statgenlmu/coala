conv_to_msms_arg <- function(feature, model) UseMethod("conv_to_msms_arg")


#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_msms_arg.default <- function(feature, model) {
  stop("Unknown feature when generating ms command")
}


#' @include simulator_class.R
#' @include simulator_ms.R
msms_class <- R6Class("Msms", inherit = ms_class,
  private = list(
    name = "msms",
    jar = NULL,
    java = NULL,
    priority = 200,
    url = "http://www.mabs.at/ewing/msms/msms3.2rc-b163.jar"
  ),
  public = list(
    initialize = function(jar, java, priority, download) {
      # Download the jar if requested by the user
      assert_that(is.logical(download) && length(download) == 1)
      if (download) {
        jar <- base::tempfile("msms_", fileext = ".jar")
        utils::download.file(private$url, jar)
      }

      # Try to automatically find a jar file and java if not given
      if (is.null(jar)) jar <- search_executable("msms.jar", "MSMS")
      if (is.null(jar)) stop("No jar file for msms found.")
      if (!file.exists(jar)) stop("msms jar (", jar, ") does not exist.")
      assert_that(is.character(jar) && length(jar) == 1)
      message("Using '", jar, "' as msms jar")
      private$jar <- jar

      if (is.null(java)) java <- search_executable(c("java", "java.exe"))
      if (is.null(java)) stop("Java not found.")
      if (!file.exists(java)) stop("Java not found.")
      assert_that(is.character(java) && length(java) == 1)
      private$java <- java

      assert_that(is.numeric(priority) && length(priority) == 1)
      private$priority <- priority

      invisible(NULL)
    },
    call = function(sample_size, n_loci, command) {
      output <- tempfile("msms")
      seed <- sample_seed(1)

      # Create the command
      cmd <- paste("-jar", private$jar,
                   paste(sample_size, n_loci, command),
                   "-seed", format(seed, scientific = FALSE))

      # Execute the command
      status <- system2(private$java, args = cmd, stdout = output)

      if (status != 0 || !file.exists(output)) stop("msms simulation failed")
      if (file.info(output)$size == 0) stop("msms output is empty")

      list(file = output, cmd = paste(private$java, cmd))
    },
    create_cmd_template = function(model) {
      cmd <- read_cache(model, "msms_cmd")
      if (is.null(cmd)) {
        cmd <- paste(vapply(model$features, conv_to_msms_arg,
                            FUN.VALUE = character(1), model),
                     collapse = "")
        cmd <- paste0("c('", cmd, "-threads 1 ')")
        cache(model, "msms_cmd", cmd)
      }
      cmd
    },
    get_info = function() {
      c(name = "msms", jar = private$jar, java = private$java)
    }
  )
)

has_msms <- function() !is.null(simulators[["msms"]])


#' Simulator: msms
#'
#' This adds the simulator 'msms' to the list of available simulators. To add
#' msms, you need to download the jar file and have Java installed on your
#' system.
#'
#' @references
#' Gregory Ewing and Joachim Hermisson.
#' MSMS: a coalescent simulation program including recombination,
#' demographic structure and selection at a single locus.
#' Bioinformatics (2010) 26 (16): 2064-2065
#' doi:10.1093/bioinformatics/btq322
#'
#' @param jar The path of the msms jar file.
#' @param java The path of the java executable on your system.
#' @param download If set to \code{TRUE}, coala will try to download
#'        the msms jar file. In that case, the \code{jar} argument
#'        is not required.
#' @inheritParams simulator_ms
#' @name simulator_msms
#' @family simulators
#' @examples
#' # Download and activate msms (requires Java)
#' \dontrun{activate_msms(download = TRUE)}
#' @export
activate_msms <- function(jar = NULL, java = NULL,
                          priority = 200, download = FALSE) {
  register_simulator(msms_class$new(jar, java,
                                    priority,
                                    download))
  reset_cache()
  invisible(NULL)
}
