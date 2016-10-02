tempfile <- function(name="unnamed") {
  base::tempfile(paste0("coala-", Sys.getpid(), "-", name, "-"))
}

require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Please install package '", pkg, "'", call. = FALSE)
  }
  invisible(TRUE)
}

sample_seed <- function(n = 1, for_ms = FALSE) {
  max_value <- ifelse(for_ms, 65536, 1e9)
  sample.int(max_value, n)
}
