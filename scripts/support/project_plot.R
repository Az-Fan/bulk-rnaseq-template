.project_plot_state <- new.env(parent = emptyenv())

load_project_plot_config <- function() {
  if (!is.null(.project_plot_state$loaded)) return(invisible(TRUE))
  if (!requireNamespace("yaml", quietly = TRUE)) stop("Package yaml is required.")

  theme_path <- Sys.getenv("THEME_CONFIG", unset = "config/theme.yml")
  registry_path <- Sys.getenv("PLOT_REGISTRY", unset = "config/plots.yml")
  .project_plot_state$theme <- yaml::read_yaml(theme_path)
  .project_plot_state$registry <- yaml::read_yaml(registry_path)
  .project_plot_state$loaded <- TRUE
  invisible(TRUE)
}

config_get <- function(x, path, default = NULL) {
  keys <- strsplit(path, "\\.", fixed = FALSE)[[1]]
  value <- x
  for (key in keys) {
    if (is.null(value[[key]])) return(default)
    value <- value[[key]]
  }
  value
}

plot_enabled <- function(path, default = TRUE) {
  load_project_plot_config()
  isTRUE(config_get(.project_plot_state$registry, path, default))
}

theme_value <- function(path, default = NULL) {
  load_project_plot_config()
  config_get(.project_plot_state$theme, path, default)
}

project_colors <- function() {
  load_project_plot_config()
  unlist(.project_plot_state$theme$colors)
}

project_theme <- function(base_size = NULL, base_family = NULL) {
  if (is.null(base_size)) base_size <- theme_value("base_size", 11)
  if (is.null(base_family)) base_family <- theme_value("font_family", "Arial")
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(color = "grey35"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold")
    )
}

save_project_plot <- function(plot, path, width = NULL, height = NULL, dpi = NULL) {
  if (is.null(width)) width <- theme_value("figure.width", 7)
  if (is.null(height)) height <- theme_value("figure.height", 5)
  if (is.null(dpi)) dpi <- theme_value("figure.dpi", 300)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = theme_value("figure.background", "white"),
    limitsize = FALSE
  )
}

write_figure_note <- function(
  figure_path,
  title,
  source,
  method,
  filters = character(),
  database = NA_character_,
  interpretation = "Exploratory visualization; biological conclusions require independent validation."
) {
  if (Sys.getenv("GENERATE_FIGURE_NOTES", unset = "0") != "1") {
    return(invisible(NULL))
  }
  note_path <- paste0(tools::file_path_sans_ext(figure_path), ".md")
  lines <- c(
    paste0("# ", title),
    "",
    paste0("- Figure: `", basename(figure_path), "`"),
    paste0("- Data source: ", source),
    paste0("- Method: ", method),
    if (!is.na(database)) paste0("- Database/gene-set source: ", database),
    if (length(filters)) paste0("- Filters: ", paste(filters, collapse = "; ")),
    paste0("- Interpretation level: ", interpretation),
    paste0("- Generated: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  )
  writeLines(lines[!is.na(lines)], note_path, useBytes = TRUE)

  registry_path <- file.path(Sys.getenv("OUTPUT_DIR", "."), "logs", "figure_registry.csv")
  row <- data.frame(
    figure = normalizePath(figure_path, winslash = "/", mustWork = FALSE),
    note = normalizePath(note_path, winslash = "/", mustWork = FALSE),
    title = title,
    source = source,
    method = method,
    database = database,
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
  dir.create(dirname(registry_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    row,
    registry_path,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(registry_path),
    append = file.exists(registry_path),
    qmethod = "double"
  )
  invisible(note_path)
}
