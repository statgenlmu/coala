executables <- new.env()

get_executable <- function(name) {
  exe <- executables[[name]]
  if (is.null(exe)) stop("No executable found for", name)
  exe
}

set_executable <- function(name, path) executables[[name]] <- path

search_executable <- function(name, file_name) {
  # See if an environment variable is given
  exe <- NULL
  exe_path <- Sys.getenv(toupper(name))
  if (exe_path != "" && file.exists(exe_path)) exe <- exe_path

  # Try to find it in the PATH folders and the Working directory
  else {
    if (Sys.info()[['sysname']] == "Windows") {
      run_path <- strsplit(Sys.getenv("PATH"), ";")[[1]]
    } else {
      run_path <- strsplit(Sys.getenv("PATH"), ":")[[1]]
    }

    candidates <- do.call(c, lapply(file_name, function(x) {
      file.path(c(getwd(), run_path), x)
    }))

    for (candidate in candidates) {
      if (file.exists(candidate)) {
        exe <- candidate
        break
      }
    }
  }

  if (is.null(exe)) return(FALSE)

  set_executable(name, exe)
  TRUE
}


has_seqgen <- function() !is.null(executables[['seqgen']])
has_msms <- function() {
  !(is.null(executables[['msms']]) || is.null(executables[['java']]))
}


.onLoad <- function(libname, pkgname) {
  search_executable("java", c("java", "java.exe"))
  search_executable("msms", "msms.jar")
  search_executable("seqgen", c("seqgen", "seq-gen",
                                "seqgen.exe", "seq-gen.exe"))
}


.onAttach <- function(libname, pkgname) {
  for (name in ls(executables)) {
    packageStartupMessage("Using ",  executables[[name]], " to provide '",
                          name, "'")
  }
}
