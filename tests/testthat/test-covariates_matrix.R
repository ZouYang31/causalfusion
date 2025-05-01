

test_that("covariates_matrix returns correct matrix", {
  df <- data.frame(
    unit = c("A", "B"),
    cov1 = c(1.1, 2.2),
    cov2 = c(3.3, 4.4)
  )

  result <- covariates_matrix(df, covariates = c("cov1", "cov2"), colname_unit = "unit")

  expect_type(result, "double")
  expect_equal(dim(result), c(2, 2))
})

test_that("covariates_matrix errors on missing arguments", {
  df <- data.frame(unit = "A", cov1 = 1)

  expect_error(covariates_matrix(), "missing")
  expect_error(covariates_matrix(df, covariates = "cov1"), "missing")
})

test_that("covariates_matrix errors if data is not a dataframe", {
  mat <- matrix(1:4, nrow = 2)
  expect_error(covariates_matrix(mat, covariates = "cov1", colname_unit = "unit"),
               "must be a data frame")
})

test_that("covariates_matrix errors if covariates are not character", {
  df <- data.frame(unit = "A", cov1 = 1)
  expect_error(covariates_matrix(df, covariates = 1, colname_unit = "unit"),
               "must be a character vector")
})

test_that("covariates_matrix errors if covariates have non-numeric columns", {
  df <- data.frame(unit = "A", cov1 = "non-numeric")
  expect_error(covariates_matrix(df, covariates = "cov1", colname_unit = "unit"),
               "must be numeric")
})

test_that("covariates_matrix errors if covariates have missing values", {
  df <- data.frame(unit = c("A", "B"), cov1 = c(1, NA))
  expect_error(covariates_matrix(df, covariates = "cov1", colname_unit = "unit"),
               "Missing values")
})



