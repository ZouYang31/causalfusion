#############################
# Data Cleaning for Data
# Preprocessing & Transformation
#############################

# -----------------------------
# Data Preprocessing Functions
# -----------------------------

#' Normalize Columns Using Min-Max Scaling
#'
#' This function applies min-max normalization to either a numeric vector or specified numeric columns
#' in a data frame. If a numeric vector is provided, it scales the values to the [0, 1] range.
#' If a data frame and column names are provided, it applies normalization only to those columns.
#'
#' @param data A numeric vector or a data frame.
#' @param columns If `data` is a data frame, a character vector specifying which columns to normalize.
#'
#' @return A normalized numeric vector (if input is numeric) or a data frame with normalized columns.
#'
#' @examples
#' # Normalize vector
#' normalize_columns(c(10, 20, 30))
#'
#' # Normalize selected columns in a data frame
#' df <- data.frame(a = 1:5, b = 6:10, c = letters[1:5])
#' normalize_columns(df, columns = c("a", "b"))
#'
#' @export
normalize_columns <- function(data, columns = NULL) {
  # Check for missing inputs
  if (missing(data) || missing(columns)) {
    stop("One or more required arguments are missing: 'data' or 'columns'.")
  }

  # Ensure 'data' is a data frame
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }

  if (!is.character(columns)) {
    stop("'covariates' must be a character vector.")
  }

  # Ensure 'columns' are all present in the data
  if (!all(columns %in% colnames(data))) {
    stop("One or more specified columns do not exist in the data.")
  }

  if (!all(sapply(data[columns], is.numeric))) {
    stop("All covariate columns must be numeric.")
  }

  # Check for NA or NaN values
  if (anyNA(data[columns]) || any(sapply(data[columns], function(col) any(is.nan(col))))) {
    stop("Missing (NA or NaN) values detected in specified columns. Please impute them before proceeding.")
  }
  # Apply min-max normalization
  for (col in columns) {
    if (is.numeric(data[[col]])) {
      min_val <- min(data[[col]], na.rm = TRUE)
      max_val <- max(data[[col]], na.rm = TRUE)
      data[[col]] <- (data[[col]] - min_val) / (max_val - min_val)
    } else {
      warning(paste("Column", col, "is not numeric and will be left unchanged."))
    }
  }

  return(data)
}
#' Reorder Rows to Prioritize a Specific Unit
#'
#' This function reorders a data frame so that rows matching a specific pattern in the `unit` column appear first.
#' This is used for placing a treated target unit at the top of the data frame
#' before applying methods like synthetic control or matching.
#'
#' @param data A data frame that contains a `unit` column.
#' @param unit_pattern A string pattern used to identify rows to bring to the top (e.g., "treated").
#'
#' @return A data frame with reordered rows, where rows matching the pattern come first.
#'
#' @examples
#' df <- data.frame(unit = c("control1", "treated", "control2"), value = c(1, 5, 2))
#' reorder_unit_first(df, "treated")
#'
#' @import dplyr
#' @export
# Function to Reorder Data, Prioritizing a Specific Unit
reorder_unit_first <- function(data, unit_pattern) {
  data %>%
    filter(grepl(unit_pattern, unit)) %>%
    bind_rows(filter(data, !grepl(unit_pattern, unit)))
}

#' Convert Data Frame to Get Outcome Matrix
#'
#' This function transforms a data frame into a numeric outcome matrix. It supports
#' both long and wide formats. If the outcome variable column is provided, the function
#' assumes the data is in long format and reshapes it to wide format. If not, it assumes
#' the data is already in wide format and simply coerces it to a matrix.
#'
#' @param data A data frame containing outcome values. Can be in wide or long format.
#' @param colname_outcome_var A character string specifying the name of the outcome variable column
#' (used in long format). If \code{NULL}, the function assumes data is already in wide format. The outcome
#' variables should not include NAs.
#' @param colname_unit A character string specifying the column name for the unit (e.g., unit ID).
#' Required if \code{colname_outcome_var} is not \code{NULL}.
#' @param colname_time A character string specifying the column name for time (e.g., year or time index).
#' Required if \code{colname_outcome_var} is not \code{NULL}.
#'
#' @return A numeric matrix where rows represent units and columns represent time periods.
#'
#' @examples
#' # Long format example
#' df_long <- data.frame(
#'   unit = rep(1:2, each = 3),
#'   time = rep(2001:2003, times = 2),
#'   outcome = c(1, 2, 3, 4, 5, 6)
#' )
#' outcome_matrix(df_long, colname_outcome_var = "outcome", colname_unit = "unit", colname_time = "time")
#'
#' # Wide format example
#' df_wide <- data.frame(
#'   t1 = c(1, 4),
#'   t2 = c(2, 5),
#'   t3 = c(3, 6)
#' )
#' outcome_matrix(df_wide, colname_outcome_var = NULL, colname_unit = NULL)
#'
#' @export
outcome_matrix <- function(data, colname_outcome_var = NULL, colname_unit = NULL, colname_time = NULL){

  # Check if the data is dataframe.
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a dataframe.")
  }

  if (is.null(colname_outcome_var)) {
    if (anyNA(data)) {
      stop("Missing values (NA) detected in the outcome column.")
    }
    # Assume already in wide format, the dataframe should only include outcome variables
    y_matrix <- as.matrix(sapply(data, as.numeric))

  } else {
    # Long format
    # Extract only outcome columns and unit, year columns.
    data <- data[, c(colname_unit, colname_time, colname_outcome_var)]
    wide_data <- reshape(data, idvar = colname_unit, timevar = colname_time, direction = "wide")
    y_matrix <- as.matrix(sapply(wide_data[, -1], as.numeric))
    if (anyNA(y_matrix)) {
      stop("Missing values (NAs) detected in the outcome column after converted into wide format.")
    }
  }
  return(y_matrix)
}


#' Extract Covariate Matrix (Z and X Matrix)
#'
#' This function aggregates the specified covariate columns by the unit column
#' and returns them as a numeric matrix, excluding the unit identifier column.
#'
#' @param data A data frame containing the covariates.
#' @param covariates A character vector specifying which covariate columns to extract.
#' @param colname_unit A string specifying the name of the unit identifier column to group by.
#'
#' @return A numeric matrix of covariates aggregated by unit.
#' @examples
#' test_data <- data.frame(
#'   unit = rep(c("A", "B"), each = 2),
#'   median_income = c(40, 50, 60, 70),
#'   median_age = c(10, 13, 15, 21),
#'   proportion = c(200, 234, 25, 32)
#' )
#' covariates <- c("median_income", "median_age", "proportion")
#' covariates_matrix(test_data, covariates, "unit")
#'
#' @export
covariates_matrix <- function(data, covariates, colname_unit) {
  if (missing(data) || missing(covariates) || missing(colname_unit)) {
    stop("One or more required arguments are missing: 'data', 'covariates', or 'colname_unit'.")
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }

  if (!is.character(covariates)) {
    stop("'covariates' must be a character vector.")
  }

  if (!all(sapply(data[covariates], is.numeric))) {
    stop("All covariate columns must be numeric.")
  }

  # Check for missing (NA) values
  if (anyNA(data[covariates]) || any(sapply(data[covariates], function(col) any(is.nan(col))))) {
    stop("Missing values detected in covariate columns. Please impute them before proceeding.")
  }

  # Aggregate using base R
  formula_str <- as.formula(paste(". ~", colname_unit))
  aggregated_df <- aggregate(formula_str, data = data[, c(colname_unit, covariates)], FUN = mean, na.rm = TRUE)
  print(aggregated_df)


  # Convert to matrix
  return(as.matrix(aggregated_df))
}


#' Compute Synthetic Weighted Averages for Selected Columns
#'
#' This function computes the synthetic group's weighted averages for the specified columns
#' within a given time-related variable (e.g., year, month, or date).
#' @param data A data frame containing a 'year' column and numeric columns of interest.
#' @param w A numeric vector of weights (should sum to 1), corresponding to control units.
#' @param covariates A character vector of column names for which to compute the weighted average.
#' @param time_value The target time (year, month, or date) for filtering the data.
#'
#' @return A one-row data frame containing the synthetic weighted averages for the selected columns.
#'
#' @examples
#' compute_synth_avg(data, w = c(0.5, 0.5), covariates = c("median_income", "median_age"), year = 2021)
#'
#' @export
compute_synth_avg <- function(data, w, covariates, time_value) {
  # Filter data to the specified year
  filtered_data <- data[data$time_value == time_value, ]

  # Assume treated unit is the first row; use only control units
  control_data <- filtered_data[-1, covariates, drop = FALSE]

  # Compute weighted averages
  synth_results <- sapply(control_data, function(col) sum(col * w))

  # Return as a one-row data frame
  return(as.data.frame(t(synth_results)))
}



# compute_synth_avg <- function(data, w) {
#
#   # Compute synthetic weighted averages
#   synth_median_income   <- sum(data[data$year == 2021, "median_income"][-1] * w)
#   synth_proportion      <- sum(data[data$year == 2021, "proportion"][-1] * w)
#   synth_median_age      <- sum(data[data$year == 2021, "median_age"][-1] * w)
#   synth_65proportion    <- sum(data[data$year == 2021, "X65p_proportion"][-1] * w)
#
#   # Store results in a named data frame for easy access
#   synth_results <- data.frame(
#     median_income  = synth_median_income,
#     proportion     = synth_proportion,
#     median_age     = synth_median_age,
#     X65p_proportion = synth_65proportion
#   )
#
#   return(synth_results)
# }


#' Compute Average Covariate Values for Control Units
#'
#' This function computes the mean of selected covariates for control units,
#' filtered by a specified time variable (e.g., year, month, or date).
#'
#' @param data A data frame containing control unit information.
#' @param covariates A character vector of covariate column names to average.
#' @param time_var A string specifying the name of the time variable column (e.g., "year", "month", or "date").
#' @param time_value A value (numeric or character) to filter the time variable by.
#'
#' @return A list containing:
#' \item{controls_avg_df}{A data frame with average covariate values for the specified time.}
#' @export
compute_controls_avg <- function(data, covariates, time_var, time_value) {
  # Filter data for the specified time condition
  data_filtered <- data[data[[time_var]] == time_value, ]

  # Compute mean for each covariate
  avg_values <- sapply(covariates, function(col) {
    mean(data_filtered[[col]], na.rm = TRUE)
  })

  # Create a data frame from the result
  controls_avg_df <- as.data.frame(as.list(avg_values))

  return(list(controls_avg_df = controls_avg_df))
}


# compute_controls_avg <- function(data){
#   # Compute averages for each covariate
#   avg_median_income   <- mean(data[data$year == 2021, "median_income"], na.rm = TRUE)
#   avg_proportion      <- mean(data[data$year == 2021, "proportion"], na.rm = TRUE)
#   avg_median_age      <- mean(data[data$year == 2021, "median_age"], na.rm = TRUE)
#   avg_X65p_proportion <- mean(data[data$year == 2021, "X65p_proportion"], na.rm = TRUE)
#
#   # Store results in a named data frame
#   controls_avg_df <- data.frame(
#     median_income   = avg_median_income,
#     proportion      = avg_proportion,
#     median_age      = avg_median_age,
#     X65p_proportion = avg_X65p_proportion
#   )
#
#   # Convert data frame to table format for display
#   controls_avg_table <- knitr::kable(controls_avg_df, format = "simple")
#
#   # Return both data frame and table
#   return(list(
#     controls_avg_df = controls_avg_df))
#
# }

# ------------------------------------------------------------------------------#
# Function: Process Data for Synthetic, Real, and Control Groups
# ------------------------------------------------------------------------------#
process_group_data <- function(df_filtered, unit_name, w) {

  # Compute average data for the filtered dataset
  avg_data <- compute_average_by_unit(df_filtered, covariate_columns)

  # Compute synthetic data using weights
  synth_data <- compute_synth_avg(avg_data, w)

  # Extract Chelsea-specific row and compute its average
  if(unit_name == "Chelsea-black"){
    chelsea_row <- df_filtered[df_filtered$unit == unit_name, ]
  }
  else {
    chelsea_row <- df_filtered[df_filtered$unit == unit_name, ]
  }
  avg_chelsea <- compute_average_by_unit(chelsea_row, covariate_columns)

  # Compute control group averages
  avg_controls <- compute_controls_avg(avg_data)

  return(list(
    avg_data = avg_data,
    synth_data = synth_data,
    avg_chelsea = avg_chelsea,
    avg_controls = avg_controls
  ))
}


#' Convert Data Frame to Numeric Matrix
#'
#' This function converts a data frame to a numeric matrix by removing non-numeric columns such as identifiers.
#'
#' @param data A data frame with at least one non-numeric column (e.g., `unit`) and the rest numeric.
#'
#' @return A numeric matrix containing only the numeric columns of the input data.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(unit = c("A", "B"), t1 = c(1, 2), t2 = c(3, 4))
#' mat <- convert_to_matrix(df)
#' }
#'
#' @import dplyr
#' @export

# Convert Data Frames to Numeric Matrices
convert_to_matrix <- function(data) {
  data %>%
    select(-unit) %>%
    as.matrix()
}

#' Generate List of B Vectors Summing to One
#'
#' This function generates a list of combinations of a budgeting vector \eqn{(b_F, b_Z, b_X)} such that their sum equals 1.
#'
#' @param step A numeric value indicating the increment step size (default is 0.01).
#' @param min_value A numeric value specifying the minimum threshold for each component (default is 0.01).
#'
#' @return A list of budgets, each represented as a list with elements `b_F`, `b_Z`, and `b_X`.
#'
#' @examples
#' B_list <- generate_b_list(step = 0.1, min_value = 0.1)
#' length(B_list)  # View number of valid combinations
#' head(B_list)
#'
#' @export
# Generate the B list
generate_b_list <- function(step = 0.01, min_value = 0.01) {
  B <- list()  # Initialize empty list
  c <- 1  # Counter

  for (b_F in seq(min_value, 1 - min_value, by = step)) {
    for (b_Z in seq(min_value, 1 - b_F, by = step)) {
      b_X <- 1 - b_F - b_Z  # Ensure sum to 1

      if (b_X >= min_value) {  # Ensure b_X meets the minimum threshold
        B[[c]] <- list(b_F = b_F, b_Z = b_Z, b_X = b_X)
        c <- c + 1
      }
    }
  }

  return(B)
}




