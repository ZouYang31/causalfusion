

test_that("outcome_matrix works for wide-format data", {
  df_wide <- data.frame(
    unit1 = c(1, 2, 3),
    unit2 = c(4, 5, 6)
  )

  result <- outcome_matrix(df_wide)

  expect_type(result, "double")
  expect_equal(dim(result), c(3, 2))
})

test_that("outcome_matrix works for long-format data", {
  df_long <- data.frame(
    unit = rep(c("A", "B"), each = 3),
    time = rep(1:3, times = 2),
    outcome = c(1, 2, 3, 4, 5, 6)
  )

  result <- outcome_matrix(df_long, colname_outcome_var = "outcome", colname_unit = "unit", colname_time = "time")

  expect_type(result, "double")
  expect_equal(dim(result), c(2, 3))
})

test_that("outcome_matrix stops on non-dataframe input", {
  non_df <- matrix(1:4, nrow = 2)

  expect_error(outcome_matrix(non_df), "must be a dataframe")
})

test_that("outcome_matrix stops if NA present in wide data", {
  df_wide <- data.frame(
    unit1 = c(1, NA, 3),
    unit2 = c(4, 5, 6)
  )

  expect_error(outcome_matrix(df_wide), "Missing values")
})

test_that("outcome_matrix stops if NA present after reshaping", {
  df_long <- data.frame(
    unit = rep(c("A", "B"), each = 3),
    time = rep(1:3, times = 2),
    outcome = c(1, 2, NA, 4, 5, 6)
  )

  expect_error(outcome_matrix(df_long, colname_outcome_var = "outcome", colname_unit = "unit", colname_time = "time"),
               "Missing values")
})
