#' @importFrom settings options_manager
CSR_OPTIONS <- options_manager(msms_path=NULL, seqgen_path=NULL)

get_msms_path <- function() CSR_OPTIONS('msms_path')
set_msms_path <- function(path) CSR_OPTIONS(msms_path=path)

get_seqgen_path <- function() CSR_OPTIONS('seqgen_path')
set_seqgen_path <- function(path) CSR_OPTIONS(seqgen_path=path)
