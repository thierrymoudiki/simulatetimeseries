
# Load required functions
source("R/univariate.R")

# Generate the time series dataset
ts_dataset <- get_data_1(diffs = FALSE)

# Create data directory if it doesn't exist
dir.create("data", showWarnings = FALSE)

# Save the dataset
usethis::use_data(ts_dataset, overwrite = TRUE)

