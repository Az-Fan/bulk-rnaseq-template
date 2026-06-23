output_dir <- Sys.getenv("PERSONALIZED_OUTPUT_DIR")
standard_results <- Sys.getenv("STANDARD_RESULTS_DIR")
standard_work <- Sys.getenv("STANDARD_WORK_DIR")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

writeLines(
  c(
    "Example personalized analysis",
    paste0("Standard results: ", standard_results),
    paste0("Standard work: ", standard_work),
    "Replace personalized/code/example_analysis.R with project-specific analysis."
  ),
  file.path(output_dir, "README.txt")
)
