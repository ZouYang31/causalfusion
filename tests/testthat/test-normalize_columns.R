test_that("normalize_columns performs min-max scaling correctly", {
  df <- data.frame(a = 1:10, b = 11:20)
  result <- normalize_columns(df, columns = c("a", "b"))
  expect_equal(result$a[1], 0)
  expect_equal(result$a[10], 1)
  expect_equal(result$b[1], 0)
  expect_equal(result$b[10], 1)
})

test_that("normalize_columns scales to expected intermediate values", {
  df <- data.frame(a = c(10, 20, 30), b = c(100, 200, 300))
  result <- normalize_columns(df, columns = c("a", "b"))

  # Expected:
  # a: (10 - 10)/(30 - 10) = 0, (20 - 10)/(30 - 10) = 0.5, (30 - 10)/(30 - 10) = 1
  # b: (100 - 100)/(300 - 100) = 0, (200 - 100)/(300 - 100) = 0.5, (300 - 100)/(300 - 100) = 1

  expect_equal(result$a, c(0, 0.5, 1))
  expect_equal(result$b, c(0, 0.5, 1))
})

test_that("normalize_columns handles uneven numeric values correctly", {
  df <- data.frame(
    x = c(3, 7, 15),
    y = c(100, 150, 250)
  )

  result <- normalize_columns(df, columns = c("x", "y"))

  # Manually compute expected normalized values:
  # For x: min = 3, max = 15 ⇒ (3−3)/(15−3)=0, (7−3)/(15−3)=4/12≈0.333, (15−3)/(15−3)=1
  # For y: min = 100, max = 250 ⇒ (100−100)/150=0, (150−100)/150=50/150≈0.333, (250−100)/150=150/150=1
  expect_equal(result$x, c(0, 1/3, 1))
  expect_equal(result$y, c(0, 1/3, 1))
})

test_that("normalize_columns returns same non-normalized columns", {
  df <- data.frame(a = 1:3, b = 4:6, c = letters[1:3])
  result <- normalize_columns(df, columns = c("a", "b"))
  expect_equal(result$c, df$c)  # Non-target column remains unchanged
})

test_that("normalize_columns throws error for missing arguments", {
  expect_error(normalize_columns(data = data.frame(a = 1), columns = NULL))
  expect_error(normalize_columns(columns = c("a")))
})

test_that("normalize_columns throws error for non-data.frame input", {
  expect_error(normalize_columns(matrix(1:9, nrow = 3), columns = c("a")))
})

test_that("normalize_columns throws error for non-character columns", {
  df <- data.frame(a = 1:3, b = 4:6)
  expect_error(normalize_columns(df, columns = 1))
})

test_that("normalize_columns throws error if column names are not in data", {
  df <- data.frame(a = 1:3, b = 4:6)
  expect_error(normalize_columns(df, columns = c("x")))
})

test_that("normalize_columns throws error if column contains NA or NaN", {
  df <- data.frame(a = c(1, NA, 3), b = c(1, 2, NaN))
  expect_error(normalize_columns(df, columns = c("a")))
  expect_error(normalize_columns(df, columns = c("b")))
})
