.onLoad <- function(libname, pkgname) {
  # scrm should always be available
  register_simulator(scrm_class$new())

  # Silently create the simulators for which a binary is available
  tryCatch(register_simulator(ms_class$new()), error = function(e) {}) #nolint
  tryCatch(register_simulator(msms_class$new()), error = function(e) {}) #nolint
  tryCatch(register_simulator(seqgen_class$new()), error = function(e) {}) #nolint
}
