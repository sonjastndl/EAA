process_pRDA <- function(formula, variables, outdir, prefix) {
  # Ensure the necessary directories exist
  ensure_directory <- function(file_path) {
    dir_path <- dirname(file_path)
    if (!file.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
      message("Directory created at: ", dir_path)
    }
  }
  
  # Define file paths
  rda_file <- paste0(outdir, "/partialRDA/pRDAobjects/", prefix, ".rds")
  plot_file <- paste0(outdir, "/partialRDA/plots/", prefix, ".png")
  rsquared_file <- paste0(outdir, "/partialRDA/Rsquared/Rsquared_", prefix, ".csv")
  anova_file <- paste0(outdir, "/partialRDA/anova/Anova_", prefix, ".csv")
  
  # Ensure directories exist
  ensure_directory(rda_file)
  ensure_directory(plot_file)
  ensure_directory(rsquared_file)
  ensure_directory(anova_file)
  
  # Check if the RDA file exists
  if (file.exists(rda_file)) {
    # Load the file
    pRDA_result <- readRDS(file = rda_file)
    message("File loaded successfully from: ", rda_file)
  } else {
    # Generate the RDA result and save it
    pRDA_result <- rda(formula, variables)
    saveRDS(pRDA_result, rda_file)
    message("File generated and saved successfully to: ", rda_file)
  }
  
  # Generate and save the plot
  png(filename = plot_file, width = 800, height = 600)
  plot(pRDA_result)
  dev.off()
  message("Plot saved successfully to: ", plot_file)
  
  # Save R-squared values
  rsquared_data <- as.data.frame(RsquareAdj(pRDA_result))
  write.csv(rsquared_data, file = rsquared_file, row.names = FALSE)
  message("R-squared values saved successfully to: ", rsquared_file)
  
  # Save ANOVA results
  anova_data <- anova(pRDA_result)
  write.csv(anova_data, file = anova_file, row.names = FALSE)
  message("ANOVA results saved successfully to: ", anova_file)
}

# Example Usage
# Replace formula_struct, Variables, outdir, and "structure" with actual values
process_pRDA(
  formula = formula_struct,
  variables = Variables,
  outdir = outdir,
  prefix = "structure"
)
