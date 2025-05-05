.onAttach <- function(libname, pkgname) {
  results_dir <- file.path("Results")
  signature_dir <- file.path(results_dir, "custom_signatures")

  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(signature_dir, showWarnings = FALSE, recursive = TRUE)
}
