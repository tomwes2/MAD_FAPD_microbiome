# Isolate MAD replication inputs vs generated outputs.
# Inputs stay in the replication kit. Generated files go to MAD_OUTPUT_DIR
# (defaults to the kit directory). Never read/write /DATA01/MAD.

.mad_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(file_arg) && nzchar(file_arg[[1]]) && file.exists(file_arg[[1]])) {
    return(dirname(normalizePath(file_arg[[1]])))
  }
  normalizePath(".", mustWork = TRUE)
}

mad_data_dir <- function() {
  env <- Sys.getenv("MAD_DATA_DIR", unset = "")
  if (nzchar(env)) return(normalizePath(env, mustWork = TRUE))
  .mad_script_dir()
}

mad_output_root <- function(data_dir = mad_data_dir()) {
  env <- Sys.getenv("MAD_OUTPUT_DIR", unset = "")
  if (nzchar(env)) {
    dir.create(env, recursive = TRUE, showWarnings = FALSE)
    return(normalizePath(env, mustWork = TRUE))
  }
  data_dir
}

mad_init_paths <- function() {
  data_dir <- mad_data_dir()
  output_root <- mad_output_root(data_dir)
  Sys.setenv(MAD_DATA_DIR = data_dir, MAD_OUTPUT_DIR = output_root)
  list(data_dir = data_dir, output_root = output_root)
}
