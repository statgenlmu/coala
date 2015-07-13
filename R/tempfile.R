tempfile <- function(name="unnamed") {
  base::tempfile(paste0("coala-", Sys.getpid(), "-", name, "-"))
}
