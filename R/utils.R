"%nin%" <- function(x, table) {
  match(x, table, nomatch = 0L) == 0L
}

check_log_file <- function(log_file) {
  if (file.exists(log_file)) {
    paste0(
      dirname(log_file),
      format(Sys.time(), "/%Y-%m-%d_%H:%M:%S_"),
      basename(log_file)
    )
  }
}

check_plot_dir <- function(plot_dir, folder_name) {

  if (!dir.exists(plot_dir)) {
    stop("The plot directory, ", plot_dir, " does not exist")
  }

  plot_folder <- paste0(plot_dir, "/", folder_name, "/")
  if (!dir.exists(plot_folder)) {
    dir.create(plot_folder)
  } else {
    stamp <- format(Sys.time(), "_%Y-%m-%d_%H:%M:%S")
    plot_folder <- paste0(plot_dir, "/", folder_name, stamp, "/")
    dir.create(plot_folder)
  }

  plot_folder

}

parse_parameters <- function(
  parameter_list,
  method = c("centwave", "align_group")
) {

  default_params <- switch(
    method,
    "centwave" = {
      param_types_centwave()
      get_centwave_defaults()
    },
    "align_group" = {
      param_types_align_group()
      c(get_obiwarp_defaults(), get_density_defaults())
    }
  )

  default_params[names(parameter_list)] <- parameter_list

  params <- list()
  params$to_optimize <- default_params

  params <- list()
  params$to_optimize <- Filter(function(x) length(x) == 2, default_params)
  params$constant <- Filter(function(x) length(x) != 2, default_params[c(quant, qual)])
  params$lists <- default_params[lists]

  params
}
